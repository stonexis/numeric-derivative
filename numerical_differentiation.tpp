#include <memory>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <utility>
#include <functional>

namespace Norms {
    enum Type : std::size_t { // Не вызывает никаких накладных расходов, поскольку на этапе компиляции преобразуется в числа
        L_1,  // 0 
        L_2,  // 1 
        L_inf, // 2 
        Count, // 3 
    };
}
/**
 * @brief Функция для генерации или измельчения массива значений заданной функции на отрезке и вычисления аналитической производной
 * @tparam T Тип данных (float, double, long double).
 * @param[out] count_nodes_out Количество узлов получившейся сетки (Выходной параметр)
 * @param x_func_derivative_rare Кортеж из трех значений старой сетки x, функции на x, производной на x (По умолчанию nullptr, nullptr, nullptr)
 * @param count_nodes_init Количество узлов для инициализации сетки, при существовании более редкой сетки ДОЛЖЕН быть равен количеству узлов в ней (По умолчанию Task_const::K)
 * @param ratio Во сколько раз необходимо измельчить сетку (По умолчанию = 1)
 * @param a Левый конец отрезка (По умолчанию Task_const::A)
 * @param b Правый конец отрезка (По умолчанию Task_const::B)
 * @return Пара массивов pair(func, antiderivative)
 * @note Каждый массив имеет размер count_nodes_out, накладные расходы по памяти - размер пары 2*размер указателя = 16 байт
 */
template <typename T>
std::tuple<T*, T*, T*> gen_grid_func_and_analytic_derivative(
                                    std::size_t& count_nodes_out,
                                    std::tuple<T*, T*, T*> x_func_derivative_rare, 
                                    const std::size_t count_nodes_init,
                                    const std::size_t ratio, 
                                    const T a, const T b
                                    ){
    auto& [x_rare, func_rare, derivative_rare] = x_func_derivative_rare; // Распаковка пары
    if (count_nodes_init < 2) throw std::invalid_argument("Invalid count_nodes_start values");                                    
    if (std::abs(b - a) < std::numeric_limits<T>::epsilon()) throw std::invalid_argument("Invalid a, b values");
    if (a > b) throw std::invalid_argument("Invalid a, b values");
    if (ratio < 1) throw std::invalid_argument("Invalid ratio");
    if (func_rare == nullptr && ratio != 1) throw std::invalid_argument("Incorrect initialization");

    T* grid_x = nullptr; // Массив сетки
    T* arr_func = nullptr; // Массив значений функции на сетке
    T* arr_derivative = nullptr; // Массив значений производной
    static auto func = [](T x) -> T {return std::sin(x); }; // Заданная функция
    static auto der = [](T x) -> T {return std::cos(x); }; // Аналитическая производная

    //Функция для измельчения сетки x
    static auto gen_grinded_grid = [](T* grid_rare, std::size_t count_nodes_grinded, std::size_t count_nodes_rare, std::size_t ratio, T step) {
        T* grid_grinded = new T[count_nodes_grinded]{};
        for (std::size_t i = 0; i < count_nodes_rare - 1; i++) {
            grid_grinded[i * ratio] = grid_rare[i]; // Старая точка
            for (std::size_t j = 1; j < ratio; j++) { // Заполнение промежуточных точек
                grid_grinded[i * ratio + j] = grid_rare[i] + j * step; 
            }
        }
        grid_grinded[count_nodes_grinded - 1] = grid_rare[count_nodes_rare - 1]; // Последняя точка
        return grid_grinded;
    };
    if (func_rare == nullptr){ // Массив другой сетки не существует, это первый запуск функции
        count_nodes_out = count_nodes_init;
        T step = std::abs(b-a) / (count_nodes_out - 1);
        grid_x = gen_uniform_grid(step, count_nodes_out, a, b); // Создаем сетку на оси x
        arr_func = new T[count_nodes_out]{};
        arr_derivative = new T[count_nodes_out]{};
        for (std::size_t i = 0; i < count_nodes_out; i++){
            arr_func[i] = func(grid_x[i]); //заданная функция
            arr_derivative[i] = der(grid_x[i]); //аналитическая производная
        }
    }
    else { // Массив старой сетки существует, нужно измельчить сетку функции
        count_nodes_out = count_nodes_init * ratio - (ratio - 1); // Проверяется на листке бумаги
        T step = std::abs(b-a) / (count_nodes_out - 1);
        grid_x = gen_grinded_grid(x_rare, count_nodes_out, count_nodes_init, ratio, step);// Измельчаем сетку x
        arr_func = new T[count_nodes_out]{};
        arr_derivative = new T[count_nodes_out]{};
        for(std::size_t i = 0; i < count_nodes_init - 1; i++) {
            arr_func[ratio * i] = func_rare[i]; // Старая точка
            arr_derivative[ratio * i] = derivative_rare[i]; //Старая точка
            for(std::size_t j = 1; j < ratio; j++) { // Промежуточные значения
                arr_func[ratio * i + j] = func(grid_x[ratio * i + j]);
                arr_derivative[ratio * i + j] = der(grid_x[ratio * i + j]);
            }
        }
        arr_func[count_nodes_out - 1] = func_rare[count_nodes_init - 1]; // Последняя точка
        arr_derivative[count_nodes_out - 1] = derivative_rare[count_nodes_init - 1];
    }
    
    return std::make_tuple(grid_x, arr_func, arr_derivative);
}
/**
 * @brief Функция для генерации равномерной сетки на отрезке [a,b]
 * @tparam T Тип данных (float, double, long double).
 * @param step - Шаг равномерной сетки. (По умолчанию Task_conts::H)
 * @param count_nodes - Количество узлов сетки. (По умолчанию Task_const::M)
 * @param a Начало отрезка. (По умолчанию Task_const::A)
 * @param b Конец отрезка. (По умолчанию Task_const::B)
 * @return T* Указатель на массив равномерной сетки.
 * @note Массив имеет размер count_nodes. (По умолчанию Task_const::M)
 */
template <typename T>
T* gen_uniform_grid(const T step, const std::size_t count_nodes, const T a, const T b) {
    if (count_nodes <= 0) throw std::invalid_argument("Invalid count_nodes values");
    if ((b - a) < std::numeric_limits<T>::epsilon()) throw std::invalid_argument("Invalid a, b values");

    T* array = new T[count_nodes]{}; 
    for (std::size_t i = 0; i < count_nodes; i++) 
        array[i] = a + step * i; // Заполняем значения, включая последний узел, равный b
    if (array[count_nodes - 1] != b)
        array[count_nodes - 1] = b;
    
    return array;
}

/**
 * @brief Функция для вычисления производной функции 
 * @tparam T Тип данных (float, double, long double).
 * @param grid_fun - Массив функции.
 * @param count_nodes - Количество узлов сетки
 * @param step - Шаг сетки
 * @return T* Указатель на массив производной функции.
 * @note Массив имеет размер count_nodes.
 */
template <typename T>
const T* gen_derivative_func(const T* grid_fun, const std::size_t count_nodes, const T step){
    if (grid_fun == nullptr) throw std::invalid_argument("grid_fun is null");
    if (count_nodes < 2) throw std::invalid_argument("Invalid count_nodes");

    T* grid_derivative = new T[count_nodes]{};
    grid_derivative[0] = (-3.0*grid_fun[0] + 4.0*grid_fun[1] - grid_fun[2]) / (2.0 * step); //Формула для самой левой точки
    for (std::size_t i = 1; i < count_nodes - 1; i++){
        grid_derivative[i] = (grid_fun[i + 1] - grid_fun[i - 1]) / (2 * step); //Формула центральных разностей
    }
    grid_derivative[count_nodes-1] = (3.0*grid_fun[count_nodes-1] - 4.0*grid_fun[count_nodes-2] + grid_fun[count_nodes-3]) / (2.0 * step);
    return grid_derivative;
}
/**
 * @brief Функция для вычисления уточненной производной функции 
 * @tparam T Тип данных (float, double, long double).
 * @param grid_derivative_more_freq Массив производной с более частым шагом 
 * @param grid_derivative_less_freq Массив производной с менее частым шагом
 * @param[out] leading_error Массив главного члена погрешности (Выходной параметр)
 * @param count_nodes_more_freq Количество узлов сетки c большей частотой (По умолчанию Task_const::M_2)
 * @param count_nodes_less_freq Количество узлов сетки с меньшей частотой (По умолчанию Task_const::M)
 * @param step_more_freq Шаг сетки с большей частотой (По умолчанию Task_const::STEP_H_2)
 * @param step_less_freq Шаг сетки с меньшей частотой (По умолчанию Task_const::H)
 * @return Пара массивов pair(runge_romberg, leading_error)
 * @note Массивы имеют размер count_nodes_less_freq. 
 */
template <typename T>
std::pair<T*, T*> gen_runge_romberg(
                const T* grid_derivative_more_freq, 
                const T* grid_derivative_less_freq,  
                const std::size_t count_nodes_more_freq,
                const std::size_t count_nodes_less_freq,
                const T step_more_freq,
                const T step_less_freq
                ){
    if (grid_derivative_more_freq == nullptr || grid_derivative_less_freq == nullptr)
        throw std::invalid_argument("grids is null");
    if (count_nodes_more_freq < 2 || count_nodes_less_freq < 2 || count_nodes_less_freq > count_nodes_more_freq)
        throw std::invalid_argument("Invalid count_nodes");
    if (std::abs(step_more_freq - step_less_freq) < std::numeric_limits<T>::epsilon())
        throw std::invalid_argument("Invalid steps");

    T* runge_romberg = new T[count_nodes_less_freq]{}; // Уточненная производная
    T* leading_error = new T[count_nodes_less_freq]{}; // Главный член погрешности
    std::size_t ratio = 2;
    for(std::size_t i = 0, j = 0; i < count_nodes_more_freq; i++){
        if (i % ratio == 0){ // j - индекс по крупной сетке, i по мелкой
            T leading_err = (grid_derivative_more_freq[i] - grid_derivative_less_freq[j]) / (ratio*ratio - 1);
            runge_romberg[j] = grid_derivative_more_freq[i] + leading_err;
            leading_error[j] = leading_err; 
            j++;
        }
    }
    return std::make_pair(runge_romberg, leading_error);
}
/**
 * @brief Функция для записи массива в файл в формате json. {"name1" : [value1, value2, ...], "name2" : [value1, value2, ...], ...}
 * @tparam T Тип данных (float, double, long double).
 * @param out Указатель на поток ввода.
 * @param length Длина массива.
 * @param name_array Имя массива.
 */
template <typename T>
void write_to_file_arr(std::ofstream &out, const T *array, std::size_t length, std::string name_array){
	if (array == nullptr || length == 0) throw std::invalid_argument("Array is null or length is zero");
 	out << "\"" << name_array << "\"" << ": ["; //форматируем файл под json
    for (std::size_t i = 0; i < length; i++) {
        out << array[i];
        if (i != length - 1) 
			out << ", ";
    }
    out << "]";
}
template <typename T>
void write_data_to_file(
                    const T *grid_M_viz,
                    const T *grid_h,
                    const T *grid_h_2,
                    const T *derivative_analytics, 
                    const T *derivative_in_h, 
                    const T *derivative_in_h_2, 
                    const T *updated_runge, 
                    const std::size_t count_h_points, 
                    const std::size_t count_h_2_points, 
                    const std::size_t count_x_points
                    ) {
	if (derivative_analytics == nullptr || derivative_in_h == nullptr || derivative_in_h_2 == nullptr || updated_runge == nullptr) 
		throw std::invalid_argument("Input derivatives cannot be null");
    if (grid_M_viz == nullptr || grid_h == nullptr || grid_h_2 == nullptr) 
		throw std::invalid_argument("Input grids cannot be null");
	if (count_h_points <= 0 || count_h_2_points <= 0 || count_x_points <= 0) 
		throw std::invalid_argument("Invalid counts points");
	
	std::ofstream out;
    out.open("data.json"); 
	if (!out.is_open()) throw std::runtime_error("Cant open file");
	else
    {
		out << '{' << '\n';
		write_to_file_arr(out, grid_M_viz, count_x_points, "grid_M_viz");
		out << ',' << '\n';
		write_to_file_arr(out, grid_h, count_h_points, "grid_h");
		out << ',' << '\n';
		write_to_file_arr(out, grid_h_2, count_h_2_points, "grid_h_2");
		out << ',' << '\n';
		write_to_file_arr(out, derivative_analytics, count_x_points, "derivative_analytics");
		out << ',' << '\n';
		write_to_file_arr(out, derivative_in_h, count_h_points, "derivative_in_h");
        out << ',' << '\n';
		write_to_file_arr(out, derivative_in_h_2, count_h_2_points, "derivative_in_h_2");
        out << ',' << '\n';
		write_to_file_arr(out, updated_runge, count_h_points, "updated_runge");
		out  << '\n' << '}';
		out.close();
    }
}
/**
 * @brief Функция для вычисления норм ошибок относительно аналитического решения.
 * @tparam T Тип данных (float, double, long double).
 * @param analytical Массив аналитической производной.
 * @param numerical Массив численной производной.
 * @param count_nodes Количество узлов сетки.
 * @return Пара массивов размера 3.
 * @note L_1_abs, L_2_abs, L_inf_abs, L_1_rel, L_2_rel, L_inf_rel
 */
template <typename T>
std::pair<T*,T*> calculate_norms(const T* analytical, const T* numerical, const std::size_t count_nodes){
    if (analytical == nullptr || numerical == nullptr)
        throw std::invalid_argument("Input arrays cannot be null.");
    if (count_nodes <= 0) 
        throw std::invalid_argument("Invalid count_nodes value.");

    T* norms_abs = new T[Norms::Count]{};
    T* norms_rel = new T[Norms::Count]{};
    T sum_abs = 0.0, sum_2_abs = 0.0, max_abs = 0.0;
    T sum_rel = 0.0, sum_2_rel = 0.0, max_rel = 0.0;
    for (std::size_t i = 0; i < count_nodes; i++) {
        T abs_error = std::abs(analytical[i] - numerical[i]);
        sum_abs += abs_error;
        sum_2_abs += abs_error * abs_error;
        max_abs = std::max(max_abs, abs_error);
        if (std::abs(analytical[i]) > std::numeric_limits<T>::epsilon()){
            T rel_error = abs_error / std::abs(analytical[i]);
            sum_rel += rel_error;
            sum_2_rel += rel_error * rel_error;
            max_rel = std::max(max_rel, rel_error);
        }
    }

    norms_abs[Norms::L_1] = sum_abs;
    norms_abs[Norms::L_2] = std::sqrt(sum_2_abs); 
    norms_abs[Norms::L_inf] = max_abs;

    norms_rel[Norms::L_1] = sum_rel; 
    norms_rel[Norms::L_2] = std::sqrt(sum_2_rel); 
    norms_rel[Norms::L_inf] = max_rel;
    
    return std::make_pair(norms_abs, norms_rel);
}
/**
 * @brief Функция для вычисления норм одного вектора.
 * @tparam T Тип данных (float, double, long double).
 * @param numerical Массив численной производной.
 * @param count_nodes Количество узлов сетки.
 * @return Массив размера 3.
 * @note L_1_abs, L_2_abs, L_inf_abs
 */
template <typename T>
T* calculate_norms(const T* numerical, const std::size_t count_nodes){
    if (numerical == nullptr) throw std::invalid_argument("Input arrays cannot be null.");
    if (count_nodes <= 0) throw std::invalid_argument("Invalid count_nodes value.");
    
    T* norms_abs = new T[Norms::Count]{};
    T sum_abs = 0.0, sum_2_abs = 0.0, max_abs = 0.0;
    for (std::size_t i = 0; i < count_nodes; i++) {
        T abs_error = std::abs(numerical[i]);
        sum_abs += abs_error;
        sum_2_abs += abs_error * abs_error;
        max_abs = std::max(max_abs, abs_error);
    }
    norms_abs[Norms::L_1] = sum_abs;
    norms_abs[Norms::L_2] = std::sqrt(sum_2_abs); 
    norms_abs[Norms::L_inf] = max_abs;
    
    return norms_abs;
}
/**
 * @brief Функция вывода значений абсолютной и относительной погрешностей в формате таблицы
 * @tparam T Тип данных (float, double, long double).
 * @param errors_h Пара массивов ошибок на сетке h
 * @param errors_h_2 Пара массивов ошибок на сетке h_2
 * @param errors_runge Пара массивов ошибок Рунге
 * @param leading_error Массив главного члена погрешности
 */
template <typename T>
void print_error_table(
                const std::pair<T*,T*>& errors_h, 
                const std::pair<T*,T*>& errors_h_2, 
                const std::pair<T*,T*>& errors_runge, 
                const T* norms_leading_error
                ){
    //Распаковка пар
    auto& [errors_h_abs, errors_h_rel] = errors_h;
    auto& [errors_h_2_abs, errors_h_2_rel] = errors_h_2;
    auto& [errors_runge_abs, errors_runge_rel] = errors_runge;
    
    // Функция для вывода строки таблицы
    static auto print_row = [](const std::string& label, const T* data) { // Первый аргумент лямбда функции - название выводимой строки, второй - данные
            std::cout << std::left // Выравнивание текста влево
                      << std::setw(18) << label // Устанавливает фиксированную ширину вывода
                      << std::scientific << std::setprecision(6) //  Устанавливает точность
                      << std::setw(15) << data[0]
                      << std::setw(15) << data[1]
                      << std::setw(15) << data[2] << "\n" ;
        };
    // Заголовок таблицы
    std::cout << std::left
              << std::setw(18) << " "
              << std::setw(15) << "L_1"
              << std::setw(15) << "L_2"
              << std::setw(15) << "L_inf" << "\n";

    std::cout << std::string(61, '-') << "\n";

    std::cout << "Errors_h" << "\n";
    print_row("Absolute error", errors_h_abs);
    print_row("Relative error", errors_h_rel);

    std::cout << std::string(61, '-') << "\n";

    std::cout << "Errors_h_2" << "\n";
    print_row("Absolute error", errors_h_2_abs);
    print_row("Relative error", errors_h_2_rel);

    std::cout << std::string(61, '-') << "\n";

    std::cout << "Errors_Runge" << "\n";
    print_row("Absolute error", errors_runge_abs);
    print_row("Relative error", errors_runge_rel);
    std::cout << std::string(61, '-') << "\n";

    print_row("Leading error", norms_leading_error);
    print_row("Absolute Errors", errors_h_abs);
    std::cout << std::string(61, '-') << "\n";
}

