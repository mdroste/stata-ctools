
#include <assert.h>
#include "ctools_hash.c"
int main(void) {
    ctools_str_hash_table ht;
    assert(ctools_str_hash_init(&ht, 16) == 0);
    const char *names[] = {"negative", "zero", "positive"};
    int codes[] = {-10, 0, 7};
    for (int i = 0; i < 3; i++) {
        assert(ctools_str_hash_insert_value(&ht, names[i], codes[i]) == 0);
        int value = 999;
        assert(ctools_str_hash_lookup(&ht, names[i], &value) == 1 && value == codes[i]);
    }
    int value = 999;
    assert(ctools_str_hash_lookup(&ht, "absent", &value) == 0 && value == 999);
    ctools_str_hash_free(&ht);
    return 0;
}
