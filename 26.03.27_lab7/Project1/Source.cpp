#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif

// ============================================================
//  Константы установки
// ============================================================
const double R3 = 12.0;   // сопротивление R3, Ом
const double TO_MA = 1000.0; // перевод А -> мА
const double TEMPERATURE = 298.0;  // T, К
const double DELTA_T = 0.5;    // погрешность T, К
const double STUDENT_T = 2.3646; // коэффициент Стьюдента, n=8, p=0.95

// ============================================================
//  Исходные данные (Таблица 1) — грубая шкала V1
//  U_eb  — напряжение эмиттер-база, В
//  U_kb  — напряжение на R3, В
// ============================================================
const std::vector<double> U_eb1 = {
    0.30, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.40, 0.41,
    0.42, 0.43, 0.44, 0.45
};

const std::vector<double> U_kb1 = {
    0.0020, 0.0041, 0.0073, 0.0125, 0.0201, 0.0285,
    0.0390, 0.0505, 0.0707, 0.0880, 0.1175, 0.1172,
    0.1760, 0.2067, 0.2320, 0.2537
};

// ============================================================
//  Исходные данные (Таблица 1') — точная шкала V1
// ============================================================
const std::vector<double> U_eb2 = {
    0.3000, 0.3108, 0.3200, 0.3297, 0.3405, 0.3497,
    0.3600, 0.3700, 0.3800, 0.3898, 0.4000, 0.4101,
    0.4202, 0.4300, 0.4401, 0.4501
};

const std::vector<double> U_kb2 = {
    0.0017, 0.0032, 0.0054, 0.0088, 0.0153, 0.0224,
    0.0333, 0.0461, 0.0627, 0.0825, 0.1057, 0.1312,
    0.1601, 0.1879, 0.2125, 0.2376
};

// ============================================================
//  Структура результата метода парных точек
// ============================================================
struct PairedPointsResult {
    std::vector<double> a_i;      // частные наклоны
    double a_mean;                // среднее tg α
    double s;                     // стандартная погрешность
    double delta_a;               // доверительный интервал tg α
    double ek;                    // e/k, К/В
    double delta_ek;              // погрешность e/k, К/В
    double ln_I0;                 // ln I0
    double I0;                    // ток насыщения, мА
};

// ============================================================
//  Вычисление коллекторного тока и его логарифма
// ============================================================
void computeCurrentAndLog(
    const std::vector<double>& U_kb_in,
    std::vector<double>& Ik,
    std::vector<double>& lnIk)
{
    int n = (int)U_kb_in.size();
    Ik.resize(n);
    lnIk.resize(n);
    for (int i = 0; i < n; ++i) {
        Ik[i] = (U_kb_in[i] / R3) * TO_MA;  // мА
        lnIk[i] = std::log(Ik[i]);
    }
}

// ============================================================
//  Метод парных точек
//  x — вектор U_eb (длина 2N)
//  y — вектор ln Ik (длина 2N)
// ============================================================
PairedPointsResult pairedPoints(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double T)
{
    int total = (int)x.size();
    int N = total / 2;   // количество пар

    PairedPointsResult res;
    res.a_i.resize(N);

    // Частные наклоны a_i = (y_{i+N} - y_i) / (x_{i+N} - x_i)
    for (int i = 0; i < N; ++i) {
        double dX = x[i + N] - x[i];
        double dY = y[i + N] - y[i];
        res.a_i[i] = dY / dX;
    }

    // Среднее
    res.a_mean = 0.0;
    for (int i = 0; i < N; ++i)
        res.a_mean += res.a_i[i];
    res.a_mean /= N;

    // Стандартная погрешность среднего  s = sqrt( sum(a_i - a_mean)^2 / (N*(N-1)) )
    double sumSq = 0.0;
    for (int i = 0; i < N; ++i) {
        double d = res.a_i[i] - res.a_mean;
        sumSq += d * d;
    }
    res.s = std::sqrt(sumSq / (N * (N - 1)));
    res.delta_a = STUDENT_T * res.s;

    // e/k и её погрешность (формула косвенных измерений)
    res.ek = T * res.a_mean;
    res.delta_ek = std::sqrt(
        (1.0 / 9.0) * (res.a_mean * res.a_mean * DELTA_T * DELTA_T
            + T * T * res.delta_a * res.delta_a)
    );

    // Ток насыщения: ln I0 = mean(ln Ik) - (e/k/T) * mean(U_eb)
    double mean_y = 0.0, mean_x = 0.0;
    for (int i = 0; i < total; ++i) {
        mean_y += y[i];
        mean_x += x[i];
    }
    mean_y /= total;
    mean_x /= total;

    res.ln_I0 = mean_y - res.a_mean * mean_x;
    res.I0 = std::exp(res.ln_I0);

    return res;
}

// ============================================================
//  Вывод таблицы 1
// ============================================================
void printTable1(
    const std::vector<double>& ueb,
    const std::vector<double>& ukb,
    const std::vector<double>& Ik,
    const std::vector<double>& lnIk)
{
    std::cout << "\n=== Таблица 1. Результаты измерений ===\n";
    std::cout << std::setw(4) << "№"
        << std::setw(10) << "U_эб, В"
        << std::setw(12) << "U_кб, В"
        << std::setw(14) << "Ik, мА"
        << std::setw(12) << "ln Ik"
        << "\n";
    std::cout << std::string(52, '-') << "\n";
    for (int i = 0; i < (int)ueb.size(); ++i) {
        std::cout << std::setw(4) << (i + 1)
            << std::fixed << std::setprecision(2)
            << std::setw(10) << ueb[i]
            << std::setprecision(4)
            << std::setw(12) << ukb[i]
            << std::setprecision(4)
            << std::setw(14) << Ik[i]
            << std::setprecision(4)
            << std::setw(12) << lnIk[i]
            << "\n";
    }
}

// ============================================================
//  Вывод таблицы 2 (метод парных точек)
// ============================================================
void printTable2(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const PairedPointsResult& res)
{
    int N = (int)res.a_i.size();
    std::cout << "\n=== Таблица 2. Метод парных точек ===\n";
    std::cout << std::setw(4) << "№"
        << std::setw(8) << "x_II"
        << std::setw(8) << "x_I"
        << std::setw(10) << "x_II-x_I"
        << std::setw(8) << "y_II"
        << std::setw(8) << "y_I"
        << std::setw(10) << "y_II-y_I"
        << std::setw(10) << "a_i"
        << std::setw(12) << "(a_i-ā)²"
        << "\n";
    std::cout << std::string(78, '-') << "\n";

    double sumA = 0.0, sumSqDev = 0.0;
    for (int i = 0; i < N; ++i) {
        double dX = x[i + N] - x[i];
        double dY = y[i + N] - y[i];
        double dev = res.a_i[i] - res.a_mean;
        sumA += res.a_i[i];
        sumSqDev += dev * dev;

        std::cout << std::fixed << std::setprecision(4)
            << std::setw(4) << (i + 1)
            << std::setw(8) << x[i + N]
            << std::setw(8) << x[i]
            << std::setw(10) << dX
            << std::setw(8) << y[i + N]
            << std::setw(8) << y[i]
            << std::setw(10) << dY
            << std::setw(10) << res.a_i[i]
            << std::setw(12) << (dev * dev)
            << "\n";
    }
    std::cout << std::string(78, '-') << "\n";
    std::cout << std::fixed << std::setprecision(4)
        << "  Σ a_i = " << sumA
        << "   Σ(a_i-ā)² = " << sumSqDev << "\n";
}

// ============================================================
//  Вывод итоговых результатов
// ============================================================
void printResults(const PairedPointsResult& res)
{
    const double EK_REF = 11604.5; // справочное значение e/k, К/В
    double error_pct = std::abs(res.ek - EK_REF) / EK_REF * 100.0;

    std::cout << "\n=== Итоговые результаты ===\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "tg α (среднее)      = " << res.a_mean << " В⁻¹\n";
    std::cout << "s (ст. погрешность) = " << res.s << " В⁻¹\n";
    std::cout << "Δ(tg α) (95%)       = " << res.delta_a << " В⁻¹\n";
    std::cout << "tg α = " << res.a_mean << " ± " << res.delta_a << " В⁻¹\n\n";

    std::cout << std::setprecision(2);
    std::cout << "e/k        = " << res.ek << " К/В\n";
    std::cout << "Δ(e/k)     = " << res.delta_ek << " К/В\n";
    std::cout << "e/k = (" << res.ek << " ± " << res.delta_ek << ") К/В\n\n";

    std::cout << std::setprecision(4);
    std::cout << "ln I0 = " << res.ln_I0 << "\n";
    std::cout << std::scientific;
    std::cout << "I0    = " << res.I0 << " мА\n\n";

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "Справочное e/k = " << EK_REF << " К/В\n";
    std::cout << "Расхождение    = " << error_pct << " %\n";
}

// ============================================================
//  Экспорт CSV для построения графиков
// ============================================================
void exportCSV(
    const std::string& filename,
    const std::vector<double>& ueb,
    const std::vector<double>& Ik,
    const std::vector<double>& lnIk,
    const PairedPointsResult& res)
{
    std::ofstream f(filename);
    f << "U_eb,Ik,lnIk,lnIk_theory\n";
    for (int i = 0; i < (int)ueb.size(); ++i) {
        double lnIk_th = res.ln_I0 + res.a_mean * ueb[i];
        f << std::fixed << std::setprecision(6)
            << ueb[i] << ","
            << Ik[i] << ","
            << lnIk[i] << ","
            << lnIk_th << "\n";
    }
    std::cout << "CSV сохранён: " << filename << "\n";
}

// ============================================================
//  Обработка одной серии: вывод всех таблиц и результатов
// ============================================================
void processSeries(
    const std::string& label,
    const std::vector<double>& ueb,
    const std::vector<double>& ukb,
    const std::string& csvFile)
{
    std::cout << "\n\n";
    std::cout << "############################################################\n";
    std::cout << "#  " << label << "\n";
    std::cout << "############################################################\n";

    // 1. Ток и логарифм
    std::vector<double> Ik, lnIk;
    computeCurrentAndLog(ukb, Ik, lnIk);

    // 2. Таблица 1 / 1'
    printTable1(ueb, ukb, Ik, lnIk);

    // 3. Метод парных точек
    PairedPointsResult res = pairedPoints(ueb, lnIk, TEMPERATURE);

    // 4. Таблица 2 / 2'
    printTable2(ueb, lnIk, res);

    // 5. Итоговые результаты
    printResults(res);

    // 6. CSV для графиков
    exportCSV(csvFile, ueb, Ik, lnIk, res);
}

// ============================================================
//  main
// ============================================================
int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    processSeries("СЕРИЯ 1  (грубая шкала V1, 2 знака)",
        U_eb1, U_kb1, "data_series1.csv");

    processSeries("СЕРИЯ 1' (точная шкала V1, 4 знака)",
        U_eb2, U_kb2, "data_series2.csv");

    return 0;
}