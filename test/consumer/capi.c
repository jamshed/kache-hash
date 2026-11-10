// A consumer of the C interface: a map from 31-mers to 32-bit counts.

#include <kache_hash/kache_capi.h>

#include <stdio.h>
#include <string.h>


int main(void)
{
    const char* const seq = "ATATCACACCCAACCTTCAAATGCCGTGCCCTAACGCCCTAATCCTGCGC";  // 20 distinct 31-mers.
    const char* seqs[1] = {seq};
    const size_t lens[1] = {strlen(seq)};

    kache_options o = kache_default_options(31, 4);
    o.expected_keys = 1000;
    o.resize_workers = 1;
    kache_table* t = NULL;
    uint32_t user = 0;
    uint64_t inserted = 0, present = 0, keys = 0;
    if(kache_new(&o, &t) != KACHE_OK || kache_register(t, &user) != KACHE_OK ||
       kache_insert(t, user, seqs, lens, 1, KACHE_ADD, 1, &inserted) != KACHE_OK ||
       kache_count(t, seqs, lens, 1, &present) != KACHE_OK || kache_size(t, &keys) != KACHE_OK)
        return 1;
    kache_free(t);

    printf("kache_capi %s, header %s: %llu 31-mers inserted, %llu present, %llu keys\n", kache_version(),
           KACHE_HASH_VERSION_STRING, (unsigned long long)inserted, (unsigned long long)present,
           (unsigned long long)keys);
    return strcmp(kache_version(), KACHE_HASH_VERSION_STRING) == 0 && present == 20 && keys == 20 ? 0 : 1;
}
