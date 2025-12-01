#include "../include/mke.h"

int save_plot_data(const char* filename) {
    FILE* f = fopen(filename, "w");
    if (!f) return 1;
    for (int i = 0; i <= kol_elem; i++) {
        double x = node[i];
        double approx = basis ? q[3*i] : q[i];
        double exact = u_func(x);
        double err = fabs(approx - exact);
        // Защита от нуля — но теперь не критична
        fprintf(f, "%.10f %.12e %.12e %.3e\n", x, approx, exact, err);
    }
    fclose(f);
    return 0;
}

int generate_plot_script(const char* script_name) {
    FILE* f = fopen(script_name, "w");
    if (!f) return 1;
    fprintf(f, "#!/bin/bash\n");
    fprintf(f, "gnuplot -p << 'EOF'\n");
    fprintf(f, "set grid\n");
    fprintf(f, "set multiplot layout 2,1 title 'МКЭ: %s базис' font ',12'\n",
            basis ? "кубический" : "линейный");
    
    // 1. Решение
    fprintf(f, "set title 'Решение'\n");
    fprintf(f, "plot 'plot.dat' u 1:2 w p pt 7 ps 1 lc 'red' t 'МКЭ', \\\n");
    fprintf(f, "     '' u 1:3 w l lw 2 lc 'blue' t 'Точное'\n");
    
    // 2. Абсолютная ошибка
    fprintf(f, "set title 'Абсолютная ошибка'\n");
    fprintf(f, "plot 'plot.dat' u 1:4 w lp pt 7 ps 0.8 lc 'purple' t '|ошибка|'\n");
    
    fprintf(f, "unset multiplot\n");
    fprintf(f, "EOF\n");
    fclose(f);
#ifndef _WIN32
    system("chmod +x plot.sh");
#endif
    return 0;
}