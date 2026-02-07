/*
 * cimport_xlsx_xml.c
 * Lightweight SAX-style XML parser for XLSX files
 */

#include "cimport_xlsx_xml.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Parser states */
typedef enum {
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

struct xlsx_xml_parser {
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
    if (!p) return NULL;

    p->callback = callback;
    p->user_data = user_data;
    p->state = STATE_TEXT;

    /* Initial text buffer */
    p->text_buf_size = 4096;
    p->text_buf = (char *)malloc(p->text_buf_size);
    if (!p->text_buf) {
        free(p);
        return NULL;
    }

    return p;
}

void xlsx_xml_parser_destroy(xlsx_xml_parser *parser)
{
    if (parser) {
        free(parser->text_buf);
        free(parser);
    }
}

static bool ensure_text_space(xlsx_xml_parser *p, size_t need)
{
    if (p->text_len + need <= p->text_buf_size) return true;

    size_t new_size = p->text_buf_size * 2;
    while (new_size < p->text_len + need) new_size *= 2;

    char *new_buf = (char *)realloc(p->text_buf, new_size);
    if (!new_buf) return false;

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
    if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES) {
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
    if (!p->callback) return true;

    p->tag_name[p->tag_name_len] = '\0';

    xlsx_xml_event event = {0};
    event.type = XLSX_XML_START_ELEMENT;
    event.tag_name_full = p->tag_name;
    event.tag_name = get_local_name(p->tag_name);
    event.attrs = p->attrs;
    event.num_attrs = p->num_attrs;
    event.is_self_closing = p->is_self_closing;

    if (!p->callback(&event, p->user_data)) {
        p->stopped = true;
        return false;
    }
    return true;
}

static bool emit_end_element(xlsx_xml_parser *p)
{
    if (!p->callback) return true;

    p->tag_name[p->tag_name_len] = '\0';

    xlsx_xml_event event = {0};
    event.type = XLSX_XML_END_ELEMENT;
    event.tag_name_full = p->tag_name;
    event.tag_name = get_local_name(p->tag_name);

    if (!p->callback(&event, p->user_data)) {
        p->stopped = true;
        return false;
    }
    return true;
}

static bool emit_text(xlsx_xml_parser *p)
{
    if (!p->callback || p->text_len == 0) return true;

    /* Skip whitespace-only text */
    bool all_whitespace = true;
    for (size_t i = 0; i < p->text_len; i++) {
        if (!isspace((unsigned char)p->text_buf[i])) {
            all_whitespace = false;
            break;
        }
    }
    if (all_whitespace) {
        p->text_len = 0;
        return true;
    }

    xlsx_xml_event event = {0};
    event.type = XLSX_XML_TEXT;
    event.text = p->text_buf;
    event.text_len = p->text_len;

    p->text_len = 0;

    if (!p->callback(&event, p->user_data)) {
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
    if (!p || p->stopped) return false;
    if (!data || len == 0) return true;

    const char *end = data + len;

    for (const char *c = data; c < end && !p->stopped; c++, p->total_pos++) {
        char ch = *c;

        switch (p->state) {
        case STATE_TEXT: {
            /* Bulk scan for next '<' using SIMD-accelerated memchr */
            const char *lt = (const char *)memchr(c, '<', (size_t)(end - c));
            if (lt == NULL) {
                /* No '<' in remaining data — copy all as text */
                size_t span = (size_t)(end - c);
                if (span > 0) {
                    if (!ensure_text_space(p, span)) {
                        set_error(p, "Out of memory");
                        return false;
                    }
                    memcpy(p->text_buf + p->text_len, c, span);
                    p->text_len += span;
                }
                p->total_pos += span - 1; /* -1 because for-loop increments */
                c = end - 1;              /* -1 because for-loop increments */
            } else if (lt == c) {
                /* '<' is the current character */
                if (!emit_text(p)) return true;
                p->state = STATE_TAG_START;
                reset_tag(p);
            } else {
                /* Copy text span before '<' */
                size_t span = (size_t)(lt - c);
                if (!ensure_text_space(p, span)) {
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
            if (ch == '/') {
                p->is_end_tag = true;
                p->state = STATE_TAG_NAME;
            } else if (ch == '!') {
                /* Could be comment <!-- or CDATA <![CDATA[ or doctype */
                p->state = STATE_DECL;
            } else if (ch == '?') {
                /* XML declaration <?...?> */
                p->state = STATE_DECL;
            } else if (isalpha((unsigned char)ch) || ch == '_') {
                p->tag_name[0] = ch;
                p->tag_name_len = 1;
                p->state = STATE_TAG_NAME;
            } else {
                set_error(p, "Invalid character after '<'");
                return false;
            }
            break;

        case STATE_TAG_NAME:
            if (isalnum((unsigned char)ch) || ch == '_' || ch == '-' || ch == ':' || ch == '.') {
                if (p->tag_name_len < XLSX_XML_MAX_TAG_LEN - 1) {
                    p->tag_name[p->tag_name_len++] = ch;
                }
            } else if (isspace((unsigned char)ch)) {
                p->state = STATE_TAG_SPACE;
            } else if (ch == '>') {
                if (p->is_end_tag) {
                    if (!emit_end_element(p)) return true;
                } else {
                    if (!emit_start_element(p)) return true;
                }
                p->state = STATE_TEXT;
            } else if (ch == '/') {
                p->is_self_closing = true;
                p->state = STATE_TAG_CLOSE;
            } else {
                set_error(p, "Invalid character in tag name");
                return false;
            }
            break;

        case STATE_TAG_SPACE:
            if (isspace((unsigned char)ch)) {
                /* skip */
            } else if (ch == '>') {
                if (p->is_end_tag) {
                    if (!emit_end_element(p)) return true;
                } else {
                    if (!emit_start_element(p)) return true;
                }
                p->state = STATE_TEXT;
            } else if (ch == '/') {
                p->is_self_closing = true;
                p->state = STATE_TAG_CLOSE;
            } else if (isalpha((unsigned char)ch) || ch == '_') {
                reset_attr(p);
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES) {
                    p->attrs[p->num_attrs].name[0] = ch;
                    p->attr_name_len = 1;
                }
                p->state = STATE_ATTR_NAME;
            } else {
                set_error(p, "Invalid character in tag");
                return false;
            }
            break;

        case STATE_ATTR_NAME:
            if (isalnum((unsigned char)ch) || ch == '_' || ch == '-' || ch == ':') {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES &&
                    p->attr_name_len < XLSX_XML_MAX_ATTR_NAME - 1) {
                    p->attrs[p->num_attrs].name[p->attr_name_len++] = ch;
                }
            } else if (ch == '=') {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES) {
                    p->attrs[p->num_attrs].name[p->attr_name_len] = '\0';
                }
                p->state = STATE_ATTR_EQ;
            } else if (isspace((unsigned char)ch)) {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES) {
                    p->attrs[p->num_attrs].name[p->attr_name_len] = '\0';
                }
                /* Wait for = */
            } else {
                set_error(p, "Invalid character in attribute name");
                return false;
            }
            break;

        case STATE_ATTR_EQ:
            if (ch == '"' || ch == '\'') {
                p->quote_char = ch;
                p->attr_value_len = 0;
                p->state = STATE_ATTR_VALUE;
            } else if (isspace((unsigned char)ch)) {
                /* skip */
            } else {
                set_error(p, "Expected quote after =");
                return false;
            }
            break;

        case STATE_ATTR_VALUE:
            if (ch == p->quote_char) {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES) {
                    p->attrs[p->num_attrs].value[p->attr_value_len] = '\0';
                    p->num_attrs++;
                }
                p->state = STATE_TAG_SPACE;
            } else {
                if (p->num_attrs < XLSX_XML_MAX_ATTRIBUTES &&
                    p->attr_value_len < XLSX_XML_MAX_ATTR_VALUE - 1) {
                    p->attrs[p->num_attrs].value[p->attr_value_len++] = ch;
                }
            }
            break;

        case STATE_TAG_CLOSE:
            if (ch == '>') {
                if (!emit_start_element(p)) return true;
                /* For self-closing, also emit end */
                if (p->is_self_closing && !p->is_end_tag) {
                    if (!emit_end_element(p)) return true;
                }
                p->state = STATE_TEXT;
            } else {
                set_error(p, "Expected '>' after '/'");
                return false;
            }
            break;

        case STATE_DECL:
            /* Skip until > for declarations/comments */
            /* Simple handling - doesn't properly parse comments with -- */
            if (ch == '>') {
                p->state = STATE_TEXT;
            }
            break;

        case STATE_COMMENT:
        case STATE_CDATA:
            /* These states aren't fully implemented - XLSX doesn't use them much */
            if (ch == '>') {
                p->state = STATE_TEXT;
            }
            break;
        }
    }

    if (is_final) {
        if (p->state != STATE_TEXT || p->text_len > 0) {
            emit_text(p);
        }
    }

    return !p->stopped;
}

const char *xlsx_xml_get_error(xlsx_xml_parser *parser)
{
    return parser ? parser->error_msg : "Invalid parser";
}

const char *xlsx_xml_get_attr(const xlsx_xml_event *event, const char *name)
{
    if (!event || !name) return NULL;

    for (int i = 0; i < event->num_attrs; i++) {
        /* Match against full name or local name (without namespace prefix) */
        const char *attr_name = event->attrs[i].name;
        const char *local_name = get_local_name(attr_name);

        if (strcmp(attr_name, name) == 0 || strcmp(local_name, name) == 0) {
            return event->attrs[i].value;
        }
    }
    return NULL;
}

size_t xlsx_xml_decode_entities(char *str, size_t len)
{
    if (!str || len == 0) return 0;

    /* Fast path: most strings contain no entities */
    if (!memchr(str, '&', len)) return len;

    char *read = str;
    char *write = str;
    char *end = str + len;

    while (read < end) {
        if (*read == '&') {
            char *amp_start = read;
            read++;

            if (read < end && *read == '#') {
                /* Numeric character reference */
                read++;
                int value = 0;
                bool is_hex = false;

                if (read < end && (*read == 'x' || *read == 'X')) {
                    is_hex = true;
                    read++;
                }

                while (read < end && *read != ';') {
                    if (is_hex) {
                        if (*read >= '0' && *read <= '9') {
                            value = value * 16 + (*read - '0');
                        } else if (*read >= 'a' && *read <= 'f') {
                            value = value * 16 + (*read - 'a' + 10);
                        } else if (*read >= 'A' && *read <= 'F') {
                            value = value * 16 + (*read - 'A' + 10);
                        }
                    } else {
                        if (*read >= '0' && *read <= '9') {
                            value = value * 10 + (*read - '0');
                        }
                    }
                    read++;
                }

                if (read < end && *read == ';') {
                    read++;
                    /* Convert to UTF-8 */
                    if (value < 0x80) {
                        *write++ = (char)value;
                    } else if (value < 0x800) {
                        *write++ = (char)(0xC0 | (value >> 6));
                        *write++ = (char)(0x80 | (value & 0x3F));
                    } else if (value < 0x10000) {
                        *write++ = (char)(0xE0 | (value >> 12));
                        *write++ = (char)(0x80 | ((value >> 6) & 0x3F));
                        *write++ = (char)(0x80 | (value & 0x3F));
                    } else {
                        *write++ = (char)(0xF0 | (value >> 18));
                        *write++ = (char)(0x80 | ((value >> 12) & 0x3F));
                        *write++ = (char)(0x80 | ((value >> 6) & 0x3F));
                        *write++ = (char)(0x80 | (value & 0x3F));
                    }
                } else {
                    /* Malformed - copy as-is */
                    while (amp_start < read) *write++ = *amp_start++;
                }
            } else {
                /* Named entity reference */
                char entity[16];
                int entity_len = 0;

                while (read < end && *read != ';' && entity_len < 15) {
                    entity[entity_len++] = *read++;
                }
                entity[entity_len] = '\0';

                if (read < end && *read == ';') {
                    read++;

                    if (strcmp(entity, "amp") == 0) {
                        *write++ = '&';
                    } else if (strcmp(entity, "lt") == 0) {
                        *write++ = '<';
                    } else if (strcmp(entity, "gt") == 0) {
                        *write++ = '>';
                    } else if (strcmp(entity, "quot") == 0) {
                        *write++ = '"';
                    } else if (strcmp(entity, "apos") == 0) {
                        *write++ = '\'';
                    } else {
                        /* Unknown entity - copy as-is */
                        while (amp_start < read) *write++ = *amp_start++;
                    }
                } else {
                    /* Malformed - copy as-is */
                    while (amp_start < read) *write++ = *amp_start++;
                }
            }
        } else {
            *write++ = *read++;
        }
    }

    *write = '\0';
    return write - str;
}

bool xlsx_xml_parse_cell_ref(const char *ref, int *col, int *row)
{
    if (!ref || !col || !row) return false;

    const char *p = ref;

    /* Skip leading $ for absolute reference */
    if (*p == '$') p++;

    /* Parse column letters */
    *col = 0;
    while (*p && isalpha((unsigned char)*p)) {
        *col = *col * 26 + (toupper((unsigned char)*p) - 'A' + 1);
        p++;
    }
    *col -= 1;  /* Convert to 0-based */

    if (*col < 0) return false;

    /* Skip $ for absolute row reference */
    if (*p == '$') p++;

    /* Parse row number */
    *row = 0;
    while (*p && isdigit((unsigned char)*p)) {
        *row = *row * 10 + (*p - '0');
        p++;
    }

    return (*row > 0);
}

int xlsx_xml_col_to_index(const char *col_str)
{
    if (!col_str) return -1;

    int index = 0;
    while (*col_str && isalpha((unsigned char)*col_str)) {
        index = index * 26 + (toupper((unsigned char)*col_str) - 'A' + 1);
        col_str++;
    }
    return index - 1;  /* 0-based */
}

void xlsx_xml_index_to_col(int index, char *buf)
{
    if (!buf) return;

    if (index < 0) {
        buf[0] = '\0';
        return;
    }

    char temp[4];
    int pos = 0;

    index++;  /* Convert to 1-based */
    while (index > 0) {
        index--;
        temp[pos++] = 'A' + (index % 26);
        index /= 26;
    }

    /* Reverse */
    for (int i = 0; i < pos; i++) {
        buf[i] = temp[pos - 1 - i];
    }
    buf[pos] = '\0';
}
