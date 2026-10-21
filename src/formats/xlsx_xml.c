#include "formats/excel_address.h"
/*
 * cimport_xlsx_xml.c
 * Lightweight SAX-style XML parser for XLSX files
 */

#include "formats/xlsx_xml.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Parser states */
typedef enum
{
    STATE_TEXT,
    STATE_TAG_START,
    STATE_TAG_NAME,
    STATE_TAG_SPACE,
    STATE_ATTR_NAME,
    STATE_ATTR_EQ,
    STATE_ATTR_VALUE,
    STATE_TAG_CLOSE,
    STATE_COMMENT,
    STATE_CDATA,
    STATE_DECL
} parser_state;

struct xlsx_xml_parser
{
    xlsx_xml_callback callback;
    void *user_data;

    parser_state state;
    bool stopped;

    /* Current element being parsed */
    char tag_name[XLSX_XML_MAX_TAG_LEN];
    int tag_name_len;
    bool is_end_tag;
    bool is_self_closing;

    /* Current attribute being parsed */
    xlsx_xml_attr attrs[XLSX_XML_MAX_ATTRIBUTES];
    int num_attrs;
    int attr_name_len;
    int attr_value_len;
    char quote_char;

    /* Text accumulator */
    char *text_buf;
    size_t text_buf_size;
    size_t text_len;

    /* Error handling */
    char error_msg[256];
    size_t total_pos;
};

/* Forward declarations */
static bool emit_start_element(xlsx_xml_parser *p);
static bool emit_end_element(xlsx_xml_parser *p);
static bool emit_text(xlsx_xml_parser *p);
static void set_error(xlsx_xml_parser *p, const char *msg);

xlsx_xml_parser *xlsx_xml_parser_create(xlsx_xml_callback callback, void *user_data)
{
    xlsx_xml_parser *p = (xlsx_xml_parser *)calloc(1, sizeof(xlsx_xml_parser));
    if (!p)
        return NULL;

    p->callback = callback;
    p->user_data = user_data;
    p->state = STATE_TEXT;

    /* Initial text buffer */
    p->text_buf_size = 4096;
    p->text_buf = (char *)malloc(p->text_buf_size);
    if (!p->text_buf)
    {
        free(p);
        return NULL;
    }

    return p;
}

void xlsx_xml_parser_destroy(xlsx_xml_parser *parser)
{
    if (parser)
    {
        free(parser->text_buf);
        free(parser);
    }
}

static bool ensure_text_space(xlsx_xml_parser *p, size_t need)
{
    if (p->text_len + need <= p->text_buf_size)
        return true;

    size_t new_size = p->text_buf_size * 2;
    while (new_size < p->text_len + need)
        new_size *= 2;

    char *new_buf = (char *)realloc(p->text_buf, new_size);
    if (!new_buf)
        return false;

    p->text_buf = new_buf;
    p->text_buf_size = new_size;
    return true;
}

static void reset_tag(xlsx_xml_parser *p)
{
    p->tag_name_len = 0;
    p->is_end_tag = false;
    p->is_self_closing = false;
    p->num_attrs = 0;
}

static void reset_attr(xlsx_xml_parser *p)
{
    if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES)
    {
        p->attrs[p->num_attrs].name[0] = '\0';
        p->attrs[p->num_attrs].value[0] = '\0';
    }
    p->attr_name_len = 0;
    p->attr_value_len = 0;
}

static const char *get_local_name(const char *full_name)
{
    const char *colon = strchr(full_name, ':');
    return colon ? (colon + 1) : full_name;
}

static bool emit_start_element(xlsx_xml_parser *p)
{
    if (!p->callback)
        return true;

    p->tag_name[p->tag_name_len] = '\0';

    xlsx_xml_event event = {0};
    event.type = XLSX_XML_START_ELEMENT;
    event.tag_name_full = p->tag_name;
    event.tag_name = get_local_name(p->tag_name);
    event.attrs = p->attrs;
    event.num_attrs = p->num_attrs;
    event.is_self_closing = p->is_self_closing;

    if (!p->callback(&event, p->user_data))
    {
        p->stopped = true;
        return false;
    }
    return true;
}

static bool emit_end_element(xlsx_xml_parser *p)
{
    if (!p->callback)
        return true;

    p->tag_name[p->tag_name_len] = '\0';

    xlsx_xml_event event = {0};
    event.type = XLSX_XML_END_ELEMENT;
    event.tag_name_full = p->tag_name;
    event.tag_name = get_local_name(p->tag_name);

    if (!p->callback(&event, p->user_data))
    {
        p->stopped = true;
        return false;
    }
    return true;
}

static bool emit_text(xlsx_xml_parser *p)
{
    if (!p->callback || p->text_len == 0)
        return true;

    /* Whitespace-only text is data too (" " cells, runs split at a space);
     * callbacks ignore text outside the elements they read. */
    xlsx_xml_event event = {0};
    event.type = XLSX_XML_TEXT;
    event.text = p->text_buf;
    event.text_len = p->text_len;

    p->text_len = 0;

    if (!p->callback(&event, p->user_data))
    {
        p->stopped = true;
        return false;
    }
    return true;
}

static void set_error(xlsx_xml_parser *p, const char *msg)
{
    snprintf(p->error_msg, sizeof(p->error_msg), "%s at position %zu", msg, p->total_pos);
}

bool xlsx_xml_parse(xlsx_xml_parser *p, const char *data, size_t len, bool is_final)
{
    if (!p || p->stopped)
        return false;
    if (!data || len == 0)
        return true;

    const char *end = data + len;

    for (const char *c = data; c < end && !p->stopped; c++, p->total_pos++)
    {
        char ch = *c;

        switch (p->state)
        {
        case STATE_TEXT:
        {
            /* Bulk scan for next '<' using SIMD-accelerated memchr */
            const char *lt = (const char *)memchr(c, '<', (size_t)(end - c));
            if (lt == NULL)
            {
                /* No '<' in remaining data — copy all as text */
                size_t span = (size_t)(end - c);
                if (span > 0)
                {
                    if (!ensure_text_space(p, span))
                    {
                        set_error(p, "Out of memory");
                        return false;
                    }
                    memcpy(p->text_buf + p->text_len, c, span);
                    p->text_len += span;
                }
                p->total_pos += span - 1; /* -1 because for-loop increments */
                c = end - 1;              /* -1 because for-loop increments */
            }
            else if (lt == c)
            {
                /* '<' is the current character */
                if (!emit_text(p))
                    return true;
                p->state = STATE_TAG_START;
                reset_tag(p);
            }
            else
            {
                /* Copy text span before '<' */
                size_t span = (size_t)(lt - c);
                if (!ensure_text_space(p, span))
                {
                    set_error(p, "Out of memory");
                    return false;
                }
                memcpy(p->text_buf + p->text_len, c, span);
                p->text_len += span;
                p->total_pos += span - 1;
                c = lt - 1; /* position at char before '<', loop will advance to '<' */
            }
            break;
        }

        case STATE_TAG_START:
            if (ch == '/')
            {
                p->is_end_tag = true;
                p->state = STATE_TAG_NAME;
            }
            else if (ch == '!')
            {
                /* Could be comment <!-- or CDATA <![CDATA[ or doctype */
                p->state = STATE_DECL;
            }
            else if (ch == '?')
            {
                /* XML declaration <?...?> */
                p->state = STATE_DECL;
            }
            else if (isalpha((unsigned char)ch) || ch == '_')
            {
                p->tag_name[0] = ch;
                p->tag_name_len = 1;
                p->state = STATE_TAG_NAME;
            }
            else
            {
                set_error(p, "Invalid character after '<'");
                return false;
            }
            break;

        case STATE_TAG_NAME:
            if (isalnum((unsigned char)ch) || ch == '_' || ch == '-' || ch == ':' || ch == '.')
            {
                if (p->tag_name_len < XLSX_XML_MAX_TAG_LEN - 1)
                {
                    p->tag_name[p->tag_name_len++] = ch;
                }
            }
            else if (isspace((unsigned char)ch))
            {
                p->state = STATE_TAG_SPACE;
            }
            else if (ch == '>')
            {
                if (p->is_end_tag)
                {
                    if (!emit_end_element(p))
                        return true;
                }
                else
                {
                    if (!emit_start_element(p))
                        return true;
                }
                p->state = STATE_TEXT;
            }
            else if (ch == '/')
            {
                p->is_self_closing = true;
                p->state = STATE_TAG_CLOSE;
            }
            else
            {
                set_error(p, "Invalid character in tag name");
                return false;
            }
            break;

        case STATE_TAG_SPACE:
            if (isspace((unsigned char)ch))
            {
                /* skip */
            }
            else if (ch == '>')
            {
                if (p->is_end_tag)
                {
                    if (!emit_end_element(p))
                        return true;
                }
                else
                {
                    if (!emit_start_element(p))
                        return true;
                }
                p->state = STATE_TEXT;
            }
            else if (ch == '/')
            {
                p->is_self_closing = true;
                p->state = STATE_TAG_CLOSE;
            }
            else if (isalpha((unsigned char)ch) || ch == '_')
            {
                reset_attr(p);
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES)
                {
                    p->attrs[p->num_attrs].name[0] = ch;
                    p->attr_name_len = 1;
                }
                p->state = STATE_ATTR_NAME;
            }
            else
            {
                set_error(p, "Invalid character in tag");
                return false;
            }
            break;

        case STATE_ATTR_NAME:
            if (isalnum((unsigned char)ch) || ch == '_' || ch == '-' || ch == ':')
            {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES &&
                    p->attr_name_len < XLSX_XML_MAX_ATTR_NAME - 1)
                {
                    p->attrs[p->num_attrs].name[p->attr_name_len++] = ch;
                }
            }
            else if (ch == '=')
            {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES)
                {
                    p->attrs[p->num_attrs].name[p->attr_name_len] = '\0';
                }
                p->state = STATE_ATTR_EQ;
            }
            else if (isspace((unsigned char)ch))
            {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES)
                {
                    p->attrs[p->num_attrs].name[p->attr_name_len] = '\0';
                }
                /* Wait for = */
            }
            else
            {
                set_error(p, "Invalid character in attribute name");
                return false;
            }
            break;

        case STATE_ATTR_EQ:
            if (ch == '"' || ch == '\'')
            {
                p->quote_char = ch;
                p->attr_value_len = 0;
                p->state = STATE_ATTR_VALUE;
            }
            else if (isspace((unsigned char)ch))
            {
                /* skip */
            }
            else
            {
                set_error(p, "Expected quote after =");
                return false;
            }
            break;

        case STATE_ATTR_VALUE:
            if (ch == p->quote_char)
            {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES)
                {
                    p->attrs[p->num_attrs].value[p->attr_value_len] = '\0';
                    p->num_attrs++;
                }
                p->state = STATE_TAG_SPACE;
            }
            else
            {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES &&
                    p->attr_value_len < XLSX_XML_MAX_ATTR_VALUE - 1)
                {
                    p->attrs[p->num_attrs].value[p->attr_value_len++] = ch;
                }
            }
            break;

        case STATE_TAG_CLOSE:
            if (ch == '>')
            {
                if (!emit_start_element(p))
                    return true;
                /* For self-closing, also emit end */
                if (p->is_self_closing && !p->is_end_tag)
                {
                    if (!emit_end_element(p))
                        return true;
                }
                p->state = STATE_TEXT;
            }
            else
            {
                set_error(p, "Expected '>' after '/'");
                return false;
            }
            break;

        case STATE_DECL:
            /* Skip until > for declarations/comments */
            /* Simple handling - doesn't properly parse comments with -- */
            if (ch == '>')
            {
                p->state = STATE_TEXT;
            }
            break;

        case STATE_COMMENT:
        case STATE_CDATA:
            /* These states aren't fully implemented - XLSX doesn't use them much */
            if (ch == '>')
            {
                p->state = STATE_TEXT;
            }
            break;
        }
    }

    if (is_final)
    {
        if (p->state != STATE_TEXT || p->text_len > 0)
        {
            emit_text(p);
        }
    }

    return !p->stopped;
}

const char *xlsx_xml_get_attr(const xlsx_xml_event *event, const char *name)
{
    if (!event || !name)
        return NULL;

    for (int i = 0; i < event->num_attrs; i++)
    {
        /* Match against full name or local name (without namespace prefix) */
        const char *attr_name = event->attrs[i].name;
        const char *local_name = get_local_name(attr_name);

        if (strcmp(attr_name, name) == 0 || strcmp(local_name, name) == 0)
        {
            return event->attrs[i].value;
        }
    }
    return NULL;
}

size_t xlsx_xml_decode_entities(char *str, size_t len)
{
    if (!str || len == 0)
        return 0;

    /* Fast path: most strings contain no entities */
    if (!memchr(str, '&', len))
        return len;

    char *read = str;
    char *write = str;
    char *end = str + len;

    while (read < end)
    {
        if (*read == '&')
        {
            char *amp_start = read;
            read++;

            if (read < end && *read == '#')
            {
                /* Numeric character reference */
                read++;
                int value = 0;
                bool is_hex = false;

                if (read < end && (*read == 'x' || *read == 'X'))
                {
                    is_hex = true;
                    read++;
                }

                while (read < end && *read != ';')
                {
                    if (is_hex)
                    {
                        if (*read >= '0' && *read <= '9')
                        {
                            value = value * 16 + (*read - '0');
                        }
                        else if (*read >= 'a' && *read <= 'f')
                        {
                            value = value * 16 + (*read - 'a' + 10);
                        }
                        else if (*read >= 'A' && *read <= 'F')
                        {
                            value = value * 16 + (*read - 'A' + 10);
                        }
                    }
                    else
                    {
                        if (*read >= '0' && *read <= '9')
                        {
                            value = value * 10 + (*read - '0');
                        }
                    }
                    if (value > 0x10FFFF)
                        value = 0x110000; /* no overflow */
                    read++;
                }

                if (read < end && *read == ';')
                {
                    read++;
                    /* Not a character (NUL, surrogate, beyond U+10FFFF): U+FFFD */
                    if (value == 0 || value > 0x10FFFF || (value >= 0xD800 && value <= 0xDFFF))
                        value = 0xFFFD;
                    /* Convert to UTF-8 */
                    if (value < 0x80)
                    {
                        *write++ = (char)value;
                    }
                    else if (value < 0x800)
                    {
                        *write++ = (char)(0xC0 | (value >> 6));
                        *write++ = (char)(0x80 | (value & 0x3F));
                    }
                    else if (value < 0x10000)
                    {
                        *write++ = (char)(0xE0 | (value >> 12));
                        *write++ = (char)(0x80 | ((value >> 6) & 0x3F));
                        *write++ = (char)(0x80 | (value & 0x3F));
                    }
                    else
                    {
                        *write++ = (char)(0xF0 | (value >> 18));
                        *write++ = (char)(0x80 | ((value >> 12) & 0x3F));
                        *write++ = (char)(0x80 | ((value >> 6) & 0x3F));
                        *write++ = (char)(0x80 | (value & 0x3F));
                    }
                }
                else
                {
                    /* Malformed - copy as-is */
                    while (amp_start < read)
                        *write++ = *amp_start++;
                }
            }
            else
            {
                /* Named entity reference */
                char entity[16];
                int entity_len = 0;

                while (read < end && *read != ';' && entity_len < 15)
                {
                    entity[entity_len++] = *read++;
                }
                entity[entity_len] = '\0';

                if (read < end && *read == ';')
                {
                    read++;

                    if (strcmp(entity, "amp") == 0)
                    {
                        *write++ = '&';
                    }
                    else if (strcmp(entity, "lt") == 0)
                    {
                        *write++ = '<';
                    }
                    else if (strcmp(entity, "gt") == 0)
                    {
                        *write++ = '>';
                    }
                    else if (strcmp(entity, "quot") == 0)
                    {
                        *write++ = '"';
                    }
                    else if (strcmp(entity, "apos") == 0)
                    {
                        *write++ = '\'';
                    }
                    else
                    {
                        /* Unknown entity - copy as-is */
                        while (amp_start < read)
                            *write++ = *amp_start++;
                    }
                }
                else
                {
                    /* Malformed - copy as-is */
                    while (amp_start < read)
                        *write++ = *amp_start++;
                }
            }
        }
        else
        {
            *write++ = *read++;
        }
    }

    *write = '\0';
    return write - str;
}

bool xlsx_xml_parse_cell_ref(const char *ref, int *col, int *row)
{
    int c, r;
    if (!col || !row || !ctools_excel_address(ref, 16384, 1048576, true, &c, &r, NULL))
        return false;
    *col = c - 1;
    *row = r;
    return true;
}

/* A non-allocating structural pass over the complete part. Keep source spans
 * rather than truncated tag-name buffers so mismatched long names also fail.
 * XML depth is bounded independently of worksheet size. */
bool xlsx_xml_validate_document(const char *data, size_t len, const char *root)
{
    if (!data || !root) return false;
    struct { const char *name; size_t len; } stack[256];
    size_t depth = 0;
    bool seen_root = false;
    const char *p = data, *end = data + len;
    if (len >= 3 && !memcmp(p, "\xef\xbb\xbf", 3)) p += 3;
    while (p < end) {
        if (*p != '<') {
            if (!depth && !isspace((unsigned char)*p)) return false;
            if (!*p) return false;
            p++;
            continue;
        }
        const char *terminator = NULL;
        size_t opener = 0, closer = 0;
        if (end - p >= 4 && !memcmp(p, "<!--", 4)) {
            terminator = "-->"; opener = 4; closer = 3;
        } else if (end - p >= 9 && !memcmp(p, "<![CDATA[", 9)) {
            if (!depth) return false;
            terminator = "]]>"; opener = 9; closer = 3;
        } else if (end - p >= 2 && !memcmp(p, "<?", 2)) {
            terminator = "?>"; opener = 2; closer = 2;
        }
        if (terminator) {
            p += opener;
            while ((size_t)(end - p) >= closer && memcmp(p, terminator, closer)) p++;
            if ((size_t)(end - p) < closer) return false;
            p += closer;
            continue;
        }
        p++;
        bool closing = p < end && *p == '/';
        if (closing) p++;
        const char *name = p;
        while (p < end && (isalnum((unsigned char)*p) || *p == '_' || *p == ':' ||
                           *p == '-' || *p == '.' || (unsigned char)*p >= 128)) p++;
        size_t n = (size_t)(p - name);
        if (!n) return false;
        bool self_closing = false;
        for (;;) {
            const char *space = p;
            while (p < end && isspace((unsigned char)*p)) p++;
            if (p == end) return false;
            if (*p == '>') { p++; break; }
            if (!closing && *p == '/' && end - p >= 2 && p[1] == '>') {
                self_closing = true; p += 2; break;
            }
            if (closing || p == space) return false;
            const char *attr = p;
            while (p < end && (isalnum((unsigned char)*p) || *p == '_' || *p == ':' ||
                               *p == '-' || *p == '.' || (unsigned char)*p >= 128)) p++;
            if (p == attr) return false;
            while (p < end && isspace((unsigned char)*p)) p++;
            if (p == end || *p++ != '=') return false;
            while (p < end && isspace((unsigned char)*p)) p++;
            if (p == end || (*p != '\'' && *p != '"')) return false;
            char quote = *p++;
            while (p < end && *p != quote) {
                if (*p == '<' || !*p) return false;
                p++;
            }
            if (p == end) return false;
            p++;
        }
        if (closing) {
            if (!depth || stack[depth - 1].len != n ||
                memcmp(stack[depth - 1].name, name, n)) return false;
            depth--;
        } else {
            if (!depth) {
                if (seen_root) return false;
                const char *local = memchr(name, ':', n);
                local = local ? local + 1 : name;
                size_t local_len = (size_t)(name + n - local);
                if (strlen(root) != local_len || memcmp(local, root, local_len)) return false;
                seen_root = true;
            }
            if (!self_closing) {
                if (depth == sizeof(stack) / sizeof(stack[0])) return false;
                stack[depth].name = name;
                stack[depth++].len = n;
            }
        }
    }
    return seen_root && depth == 0;
}
