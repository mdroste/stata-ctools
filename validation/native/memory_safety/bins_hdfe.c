/* Compare batch projections with an independent dummy-variable least-squares solve. */
#include "alloc.h"
#include "stplugin.h"
#include <math.h>
#include "cbinscatter/cbinscatter_resid.c"
#include "ctools_ols.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
static ST_retcode message(char *s)
{
    (void)s;
    return 0;
}
enum
{
    OBS = 48,
    COLS = 5,
    DUMMIES = 6
};

static void reference(const double *input, const int *fe, const double *weights, double *out)
{
    double gram[DUMMIES][DUMMIES + 1] = {{0}};
    for (int i = 0; i < OBS; i++)
    {
        double row[DUMMIES] = {
            1, fe[i] == 2, fe[i] == 3, fe[OBS + i] == 2, fe[OBS + i] == 3, fe[OBS + i] == 4};
        double w = weights ? weights[i] : 1;
        for (int j = 0; j < DUMMIES; j++)
        {
            for (int k = 0; k < DUMMIES; k++)
                gram[j][k] += w * row[j] * row[k];
            gram[j][DUMMIES] += w * row[j] * input[i];
        }
    }
    for (int j = 0; j < DUMMIES; j++)
    {
        double pivot = gram[j][j];
        assert(fabs(pivot) > 1e-6);
        for (int k = j; k <= DUMMIES; k++)
            gram[j][k] /= pivot;
        for (int l = 0; l < DUMMIES; l++)
        {
            if (l == j)
                continue;
            double factor = gram[l][j];
            for (int k = j; k <= DUMMIES; k++)
                gram[l][k] -= factor * gram[j][k];
        }
    }
    for (int i = 0; i < OBS; i++)
    {
        out[i] = input[i] - gram[0][DUMMIES];
        if (fe[i] > 1)
            out[i] -= gram[fe[i] - 1][DUMMIES];
        if (fe[OBS + i] > 1)
            out[i] -= gram[1 + fe[OBS + i]][DUMMIES];
    }
}

int main(int argc, char **argv)
{
    host.spoutsml = host.spouterr = message;
    int failure = argc > 1 ? atoi(argv[1]) : 0;
    int weighted = argc > 2 ? atoi(argv[2]) : 0;
    int first = argc > 3 ? atoi(argv[3]) : 0;
    int maxiter = argc > 4 ? atoi(argv[4]) : 10000;
    int fe[2 * OBS];
    double weights[OBS], values[COLS][OBS], expected[COLS][OBS], *vars[COLS];
    for (int i = 0; i < OBS; i++)
    {
        fe[i] = 1 + (i / 16) % 3;
        fe[OBS + i] = 1 + (i / 12) % 4;
        weights[i] = 1 + i % 5;
        values[0][i] = 0;         /* The old stopping rule immediately accepts this column. */
        values[1][i] = i % 7 < 3; /* Bin indicator with nonzero FE group means. */
        values[2][i] = sin(0.7 * i) + i * 0.05;
        values[3][i] = cos(1.1 * i);
        values[4][i] = 2 * fe[i] - 3 * fe[OBS + i]; /* Fully absorbed column. */
    }
    const double *w = weighted ? weights : NULL;
    /* Also test a nonzero first column already orthogonal to every FE. */
    if (first == 1)
    {
        reference(values[3], fe, w, expected[3]);
        memcpy(values[0], expected[3], sizeof(values[0]));
    }
    for (int v = 0; v < COLS; v++)
    {
        reference(values[v], fe, w, expected[v]);
        vars[v] = values[(v + (first == 2)) % COLS];
    }
    fail_at = failure;
    int dropped = -1;
    int rc =
        hdfe_residualize_batch(vars, COLS, fe, OBS, 2, w, weighted, maxiter, 1e-12, 0, &dropped);
    double error = 0, group_mean = 0;
    for (int v = 0; v < COLS; v++)
    {
        for (int i = 0; i < OBS; i++)
            error = fmax(error, fabs(values[v][i] - expected[v][i]));
        for (int g = 0; g < 2; g++)
        {
            for (int level = 1; level <= 4; level++)
            {
                double sum = 0, count = 0;
                for (int i = 0; i < OBS; i++)
                {
                    if (fe[g * OBS + i] != level)
                        continue;
                    double wi = w ? w[i] : 1;
                    sum += wi * values[v][i];
                    count += wi;
                }
                if (count)
                    group_mean = fmax(group_mean, fabs(sum / count));
            }
        }
    }
    printf("failure=%d rc=%d allocations=%d alive=%d error=%.9g group_mean=%.9g\n", failure, rc,
           calls, alive, error, group_mean);
    leftovers();
    fflush(stdout);
    assert(alive == 0);
    if (failure)
        assert(rc == CBINSCATTER_ERR_MEMORY);
    else if (maxiter <= 1)
        assert(rc == CBINSCATTER_ERR_HDFE);
    else
    {
        assert(rc == 0 && dropped == 0);
        assert(error < 1e-9 && group_mean < 1e-10);
    }
    return 0;
}
