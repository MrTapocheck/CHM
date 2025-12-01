#include "../include/mke.h"

int save_plot_data(const char* filename) {
    FILE* f = fopen(filename, "w");
    if (!f) return 1;
    for (int i = 0; i <= kol_elem; i++) {
        double x = node[i];
        double approx = basis ? q[3*i] : q[i];
        double exact = u_func(x);
        double err = fabs(approx - exact);
        double safe_err = (err < 1e-30) ? 1e-30 : err;
        // fprintf(f, "%.10f %.12e %.12e %.12e\n", x, approx, exact, safe_err);

        fprintf(f, "%.10f\t%.12e\t%.12e\t%.12e\n", x, approx, exact, safe_err);  // ← \t
    }
    fclose(f);
    return 0;
}

int generate_plot_script(const char* script_name) {
    FILE* f = fopen(script_name, "w");
    if (!f) return 1;
    fprintf(f, "#!/bin/bash\n");
    fprintf(f, "gnuplot -p << 'EOF'\n");
    fprintf(f, "set datafile separator \"\t\"\n");  // ← настоящий таб
    // fprintf(f, "set datafile separator \"\\t\"\n");  // ← табы
    fprintf(f, "set grid\n");
    fprintf(f, "set multiplot layout 3,1 title 'МКЭ: %s базис' font ',12'\n",
            basis ? "кубический" : "линейный");
    fprintf(f, "set title 'Решение'\n");
    fprintf(f, "plot 'plot.dat' u 1:2 w p pt 7 ps 1 lc 'red' t 'МКЭ', \\\n");
    fprintf(f, "     '' u 1:3 w l lw 2 lc 'blue' t 'Точное'\n");
    fprintf(f, "set title 'Абсолютная ошибка'\n");
    fprintf(f, "unset logscale y\n");
    fprintf(f, "plot 'plot.dat' u 1:4 w lp pt 7 ps 0.7 lc 'purple' t '|ош|'\n");
    fprintf(f, "set title 'log_{10}(|ош|)'\n");
    fprintf(f, "set logscale y\n");
    fprintf(f, "plot 'plot.dat' u 1:(log10(column(4))) w lp pt 7 ps 0.7 lc 'green' t 'log_{10}(|ош|)'\n");
    fprintf(f, "unset multiplot\n");
    fprintf(f, "EOF\n");
    fclose(f);
#ifndef _WIN32
    system("chmod +x plot.sh");
#endif
    return 0;
}