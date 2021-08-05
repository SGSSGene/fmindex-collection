/*
 * Run Queries
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>

#include <sdsl/suffix_arrays.hpp>

#include "interface.h"

/* only for getTime() */
#include <sys/resource.h>
#include <sys/time.h>

#define COUNT ('C')
#define LOCATE ('L')
#define EXTRACT ('E')
#define DISPLAY ('D')
#define VERBOSE ('V')

using namespace sdsl;
using namespace std;

/* local headers */
void do_count(const CSA_TYPE &);
void do_locate(const CSA_TYPE &);
void pfile_info(ulong * length, ulong * numpatt);
double getTime(void);
void usage(char * progname);

static int Verbose = 0;
static ulong Index_size, Text_length;
static double Load_time;

/*
 * Temporary usage: run_queries <index file> <type> [length] [V]
 */
int main(int argc, char * argv[])
{
    char * filename;
    char querytype;

    if (argc < 2)
    {
        usage(argv[0]);
        exit(1);
    }

    filename = argv[1];
    querytype = *argv[2];
    string index_file = string(argv[1]) + "." + string(SUF);
    bool mapped = false;
#ifdef USE_HP
    uint64_t index_size = util::file_size(index_file);
    // allocate hugepages, add 10 MiB + 3% overhead
    uint64_t alloc_size = (uint64_t)(index_size * 1.03 + (1ULL << 20) * 10);
    try
    {
        memory_manager::use_hugepages(alloc_size);
        mapped = true;
    }
    catch (...)
    {
        std::cout << "Unable to allocate " << alloc_size << " bytes ";
        std::cout << "of hugepage space on your system" << std::endl;
    }
#endif

    CSA_TYPE csa;
    fprintf(stderr, "# File = %s\n", (string(filename) + "." + string(SUF)).c_str());
    fprintf(stderr, "# program = %s\n", string(SUF).c_str());
    Load_time = getTime();
    load_from_file(csa, index_file);
    Load_time = getTime() - Load_time;
    fprintf(stderr, "# Load_index_time_in_sec = %.2f\n", Load_time);
    std::cerr << "# text_size = " << csa.size() - 1 << std::endl;

    Index_size = size_in_bytes(csa);
    Text_length = csa.size() - 1; // -1 since we added a sentinel character
    /*	Index_size /=1024; */
    fprintf(stderr, "# Index_size_in_bytes = %lu\n", Index_size);

    fprintf(stderr, "# hugepages = %i\n", (int)mapped);
    bool use_sse = false;
#ifdef __SSE4_2__
    use_sse = true;
#endif
    fprintf(stderr, "# sse = %i\n", (int)use_sse);
    bool use_popcount_tl = false;
#ifdef POPCOUNT_TL
    use_popcount_tl = true;
#endif
    fprintf(stderr, "# popcount_tl = %i\n", (int)use_popcount_tl);

    switch (querytype)
    {
        case COUNT:
            if (argc > 3)
                if (*argv[3] == VERBOSE)
                {
                    Verbose = 1;
                    fprintf(stdout, "%c", COUNT);
                }
            do_count(csa);
            break;
        case LOCATE:
            if (argc > 3)
                if (*argv[3] == VERBOSE)
                {
                    Verbose = 1;
                    fprintf(stdout, "%c", LOCATE);
                }
            do_locate(csa);
            break;
        default: fprintf(stderr, "Unknow option: main ru\n"); exit(1);
    }
    return 0;
}

void do_count(const CSA_TYPE & csa)
{
    ulong numocc, length, tot_numocc = 0, numpatt, res_patt;
    double time, tot_time = 0;
    uchar * pattern;

    pfile_info(&length, &numpatt);
    res_patt = numpatt;

    pattern = (uchar *)malloc(sizeof(uchar) * (length));
    if (pattern == NULL)
    {
        fprintf(stderr, "Error: cannot allocate\n");
        exit(1);
    }

    while (res_patt)
    {

        if (fread(pattern, sizeof(*pattern), length, stdin) != length)
        {
            fprintf(stderr, "Error: cannot read patterns file\n");
            exit(1);
        }

        /* Count */
        time = getTime();
        numocc = sdsl::count(csa, pattern, pattern + length);

        if (Verbose)
        {
            fwrite(&length, sizeof(length), 1, stdout);
            fwrite(pattern, sizeof(*pattern), length, stdout);
            fwrite(&numocc, sizeof(numocc), 1, stdout);
        }
        tot_time += (getTime() - time);
        tot_numocc += numocc;
        res_patt--;
    }

    fprintf(stderr, "# Total_Num_occs_found = %lu\n", tot_numocc);
    fprintf(stderr, "# Count_time_in_milli_sec = %.4f\n", tot_time * 1000);
    fprintf(stderr, "# Count_time/Pattern_chars = %.4f\n", (tot_time * 1000) / (length * numpatt));
    fprintf(stderr, "# Count_time/Num_patterns = %.4f\n\n", (tot_time * 1000) / numpatt);
    fprintf(stderr,
            "# (Load_time+Count_time)/Pattern_chars = %.4f\n",
            ((Load_time + tot_time) * 1000) / (length * numpatt));
    fprintf(stderr, "# (Load_time+Count_time)/Num_patterns = %.4f\n\n", ((Load_time + tot_time) * 1000) / numpatt);

    free(pattern);
}

void do_locate(const CSA_TYPE & csa)
{
    ulong numocc, length;
    ulong tot_numocc = 0, numpatt;
    double time, tot_time = 0;
    uchar * pattern;

    pfile_info(&length, &numpatt);

    pattern = (uchar *)malloc(sizeof(uchar) * (length));
    if (pattern == NULL)
    {
        fprintf(stderr, "Error: cannot allocate\n");
        exit(1);
    }

    while (numpatt)
    {

        if (fread(pattern, sizeof(*pattern), length, stdin) != length)
        {
            fprintf(stderr, "Error: cannot read patterns file\n");
            perror("run_queries");
            exit(1);
        }
        // Locate
        time = getTime();
        auto occ = sdsl::locate(csa, (char *)pattern, (char *)pattern + length);
        numocc = occ.size();
        tot_time += (getTime() - time);

        tot_numocc += numocc;
        numpatt--;

        if (Verbose)
        {
            fwrite(&length, sizeof(length), 1, stdout);
            fwrite(pattern, sizeof(*pattern), length, stdout);
            fwrite(&numocc, sizeof(numocc), 1, stdout);
        }
    }

    fprintf(stderr, "# Total_Num_occs_found = %lu\n", tot_numocc);
    fprintf(stderr, "# Locate_time_in_secs = %.2f\n", tot_time);
    fprintf(stderr, "# Locate_time/Num_occs = %.4f\n\n", (tot_time * 1000) / tot_numocc);
    fprintf(stderr, "# (Load_time+Locate_time)/Num_occs = %.4f\n\n", ((tot_time + Load_time) * 1000) / tot_numocc);

    free(pattern);
}

/* Open patterns file and read header */
void pfile_info(ulong * length, ulong * numpatt)
{
    int error;
    uchar c;
    uchar origfilename[257];

    error = fscanf(stdin, "# number=%lu length=%lu file=%s forbidden=", numpatt, length, origfilename);
    if (error != 3)
    {
        fprintf(stderr, "Error: Patterns file header not correct\n");
        perror("run_queries");
        exit(1);
    }

    fprintf(stderr, "# pat_cnt = %lu\n", *numpatt);
    fprintf(stderr, "# pat_length = %lu\n", *length);
    fprintf(stderr, "# forbidden_chars = ");

    while ((c = fgetc(stdin)) != 0)
    {
        if (c == '\n') break;
        fprintf(stderr, "%d", c);
    }

    fprintf(stderr, "\n");
}

double getTime(void)
{

    double usertime, systime;
    struct rusage usage;

    getrusage(RUSAGE_SELF, &usage);

    usertime = (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1000000.0;

    systime = (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1000000.0;

    return (usertime + systime);
}

void usage(char * progname)
{
    fprintf(stderr, "\nThe program loads <index> and then executes over it the\n");
    fprintf(stderr, "queries it receives from the standard input. The standard\n");
    fprintf(stderr, "input comes in the format of the files written by \n");
    fprintf(stderr, "genpatterns or genintervals.\n");
    fprintf(stderr, "%s reports on the standard error time statistics\n", progname);
    fprintf(stderr, "regarding to running the queries.\n\n");
    fprintf(stderr, "Usage:  %s <index> <type> [length] [V]\n", progname);
    fprintf(stderr, "\n\t<type>   denotes the type of queries:\n");
    fprintf(stderr, "\t         %c counting queries;\n", COUNT);
    fprintf(stderr, "\t         %c locating queries;\n", LOCATE);
    fprintf(stderr, "\t         %c displaying queries;\n", DISPLAY);
    fprintf(stderr, "\t         %c extracting queries.\n\n", EXTRACT);
    fprintf(stderr, "\n\t[length] must be provided in case of displaying queries (D)\n");
    fprintf(stderr, "\t         and denotes the number of characters to display\n");
    fprintf(stderr, "\t         before and after each pattern occurrence.\n");
    fprintf(stderr, "\n\t[V]      with this options it reports on the standard output\n");
    fprintf(stderr, "\t         the results of the queries. The results file should be\n");
    fprintf(stderr, "\t         compared with trusted one by compare program.\n\n");
}
