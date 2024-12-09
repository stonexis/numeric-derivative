#include <memory>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <utility>
#include <functional>
#include <map>

/**
 * @brief Функция для генерации массива значений заданной функции на отрезке
 * @tparam T Тип данных (float, double, long double).
 * @param mode True - задание функции, False - задание производной (аналитически)
 * @param array_x Массив, для которого генерируются значения функции (Массив случайных точек).
 * @param length Длина массива точек.
 * @return T* Указатель на массив значений функции.
 * @note Имеет размер length.
 */
template <typename T>
T* gen_func_arr(bool mode, const T *array_x, const std::size_t length) {
    if (array_x == nullptr) throw std::invalid_argument("array_x is null");
    if (length <= 0) throw std::invalid_argument("Invalid length values");

    T* arr_func = new T[length]{};
    for (std::size_t i = 0; i < length; i++){
        if (mode == true)
            arr_func[i] = sin(array_x[i]); //заданная функция
        else
            arr_func[i] = cos(array_x[i]); //аналитическая производная
    }
    return arr_func;
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
const T* gen_uniform_grid(const T step, const std::size_t count_nodes, const T a, const T b) {
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
    grid_derivative[0] = (-3*grid_fun[0] + 4*grid_fun[1] - grid_fun[2]) / (2 * step); //Формула для самой левой точки
    for (std::size_t i = 1; i < count_nodes - 1; i++){
        grid_derivative[i] = (grid_fun[i + 1] - grid_fun[i - 1]) / (2 * step); //Формула центральных разностей
    }
    grid_derivative[count_nodes-1] = (3*grid_fun[count_nodes-1] - 4*grid_fun[count_nodes-2] + grid_fun[count_nodes-3]) / (2 * step);
    return grid_derivative;
}
/**
 * @brief Функция для вычисления уточненной производной функции 
 * @tparam T Тип данных (float, double, long double).
 * @param grid_derivative_more_freq Массив производной с более частым шагом 
 * @param grid_derivative_less_freq Массив производной с менее частым шагом
 * @param[out] leading_error Массив главного члена погрешности (Выходной параметр)
 * @param count_nodes_great_freq Количество узлов сетки c большей частотой (По умолчанию Task_const::M_2)
 * @param count_nodes_small_freq Количество узлов сетки с меньшей частотой (По умолчанию Task_const::M)
 * @param step_greater_freq Шаг сетки с большей частотой (По умолчанию Task_const::STEP_H_2)
 * @param step_smaller_freq Шаг сетки с меньшей частотой (По умолчанию Task_const::H)
 * @return T* Указатель на массив уточненной производной функции.
 * @note Массив имеет размер count_nodes_great_freq. (По умолчанию Task_const::M_2)
 */
template <typename T>
const T* gen_runge_romberg(
                const T* grid_derivative_more_freq, 
                const T* grid_derivative_less_freq,
                T **leading_error,  
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

    T* runge_romberg = new T[count_nodes_more_freq]{}; // Уточненная производная
    *leading_error = new T[count_nodes_less_freq]{}; // Главный член погрешности
    std::size_t ratio = 2;
    for(std::size_t i = 0, j = 0; i < count_nodes_more_freq; i++){

        if (i % ratio == 0){
            T leading_err = (grid_derivative_more_freq[i] - grid_derivative_less_freq[j]) / (ratio*ratio - 1);
            runge_romberg[i] = grid_derivative_more_freq[i] + leading_err;
            (*leading_error)[j] = leading_err; 
            j++;
        }
        else
            runge_romberg[i] = grid_derivative_more_freq[i];
    }
    return runge_romberg;
}

/**
 * @brief Функция для генерации более мелкой, равномерной сетки, внутри существующей
 * @tparam T Тип данных (float, double, long double).
 * @param arr_old Интервал, внутри которого строится мелкая сетка
 * @param length_old Размер массива, внутри которого строится сетка
 * @param length_new Размер массива новой сетки 
 * @param count_new_nodes Количество добавляемых точек (По умолчанию Task_const::M_viz)
 * @return T* Указатель на массив новой равномерной сетки.
 * @note Массив имеет размер length_new.
 */
template <typename T>
T* gen_uniform_grid_in_exist(const T* arr_old, const std::size_t length_old, std::size_t& length_new, const std::size_t count_new_nodes) {
    if (arr_old == nullptr) throw std::invalid_argument("array_old is null");
    if (length_old < 2) throw std::invalid_argument("length_old must be at least 2");
    if (count_new_nodes < 2) throw std::invalid_argument("count_new_nodes must be at least 2");

    T step = (abs(arr_old[length_old-1] - arr_old[0]))/(count_new_nodes - 1);
    length_new = length_old + count_new_nodes - 2; // Учёт старых узлов и новых точек, начало и конец отрезка содержится в обоих сетках, поэтому отнимаем 2
    T* arr_new = new T[length_new]{};
    std::copy(arr_old, arr_old + length_old, arr_new); // Копируем старую сетку
    // Добавляем новые точки равномерно, с учётом step
    T a = arr_old[0];
    std::size_t insert_index = length_old; // Индекс для вставки новых точек
    for (std::size_t i = 1; i < count_new_nodes-1; i++) {
        T value = a + step * i;
        auto position = std::lower_bound(arr_new, arr_new + insert_index, value); // Найти позицию для вставки
        std::rotate(position, arr_new + insert_index, arr_new + insert_index + 1); // Сдвинуть элементы
        *position = value; // Вставить новую точку
        insert_index++;
    }
    return arr_new;
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
		write_to_file_arr(out, updated_runge, count_h_2_points, "updated_runge");
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
 * @return Массив размера 3*2.
 * @note L_1_abs, L_2_abs, L_inf_abs, L_1_rel, L_2_rel, L_inf_rel
 */
template <typename T>
T* calculate_norms(const T* analytical, const T* numerical, const std::size_t count_nodes){
    if (analytical == nullptr || numerical == nullptr)
        throw std::invalid_argument("Input arrays cannot be null.");
    if (count_nodes <= 0) 
        throw std::invalid_argument("Invalid count_nodes value.");
    std::size_t count_norms = 3; //Количество норм для вычисления (L_1,L_2,L_inf)
    std::size_t count_errors = 2; //Количество невязок для вычисления (absolute,relative)
    T* norms = new T[count_norms * count_errors]{};
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
    T L_1_abs = sum_abs, L_2_abs = std::sqrt(sum_2_abs), L_inf_abs = max_abs;
    T L_1_rel = sum_rel, L_2_rel = std::sqrt(sum_2_rel), L_inf_rel = max_rel;
    norms[0] = L_1_abs; norms[1] = L_2_abs; norms[2] = L_inf_abs;
    norms[3] = L_1_rel; norms[4] = L_2_rel; norms[5] = L_inf_rel;
    
    return norms;
}
/**
 * @brief Функция для вычисления норм.
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
    std::size_t count_norms = 3; //Количество норм для вычисления (L_1,L_2,L_inf)
    T* norms = new T[count_norms]{};
    T sum_abs = 0.0, sum_2_abs = 0.0, max_abs = 0.0;
    for (std::size_t i = 0; i < count_nodes; i++) {
        T abs_error = std::abs(numerical[i]);
        sum_abs += abs_error;
        sum_2_abs += abs_error * abs_error;
        max_abs = std::max(max_abs, abs_error);
    }
    T L_1_abs = sum_abs, L_2_abs = std::sqrt(sum_2_abs), L_inf_abs = max_abs;
    norms[0] = L_1_abs; norms[1] = L_2_abs; norms[2] = L_inf_abs;
    
    return norms;
}

/**
 * @brief Функция вывода значений абсолютной и относительной погрешностей в формате таблицы
 * @tparam T Тип данных (float, double, long double).
 * @param errors_h Массив ошибок на сетке h
 * @param errors_h_2 Массив ошибок на сетке h_2
 * @param errors_runge Массив ошибок Рунге
 * @param leading_error Массив главного члена погрешности
 * @param count_nodes_h Количество узлов в сетке h
 */
template <typename T>
void print_error_table(const T* errors_h, const T* errors_h_2, const T* errors_runge, const T* leading_error, const std::size_t count_nodes_h){
    if (errors_h == nullptr || errors_h_2 == nullptr || errors_runge == nullptr)
        throw std::invalid_argument("Input arrays cannot be null");
    // Функция для вывода 
    // Тип функции явно указан через std::function, можно заметить на auto
    std::function<void(const std::string&, const std::string&, const T*)> print_row = 
        [](const std::string& label_first, const std::string& label_last, const T* data) { // Первый аргумент лямбда функции - названии выводимой строки, второй - данные
            std::cout << std::left // Выравнивание текста влево
                      << std::setw(18) << label_first // Устанавливает фиксированную ширину вывода
                      << std::scientific << std::setprecision(6) //  Устанавливает точность
                      << std::setw(15) << data[0]
                      << std::setw(15) << data[1]
                      << std::setw(15) << data[2] << "\n" ;
            std::cout << std::left // Выравнивание текста влево
                      << std::setw(18) << label_last // Устанавливает фиксированную ширину вывода
                      << std::scientific << std::setprecision(6) //  Устанавливает точность
                      << std::setw(15) << data[3]
                      << std::setw(15) << data[4]
                      << std::setw(15) << data[5] << "\n" ;
        };
    // Заголовок таблицы
    std::cout << std::left
              << std::setw(18) << " "
              << std::setw(15) << "L_1"
              << std::setw(15) << "L_2"
              << std::setw(15) << "L_inf" << "\n";

    std::cout << std::string(61, '-') << "\n";
    std::cout << "Errors_h" << "\n";
    print_row("Absolute error", "Relative error", errors_h);

    std::cout << std::string(61, '-') << "\n";
    std::cout << "Errors_h_2" << "\n";
    print_row("Absolute error", "Relative error", errors_h_2);

    std::cout << std::string(61, '-') << "\n";
    std::cout << "Errors_Runge" << "\n";
    print_row("Absolute error", "Relative error", errors_runge);
    std::cout << std::string(61, '-') << "\n";
    compare_errors(leading_error, errors_h, count_nodes_h);
}

/**
 * @brief Функция для сравнения норм главного члена погрешности с абсолютными ошибками
 * @tparam T Тип данных (float, double, long double).
 * @param leading_error Массив главного члена погрешности.
 * @param absolute_errors_h Нормы абсолютной погрешности на сетке с шагом h.
 * @param count_nodes Количество узлов сетки с шагом h.
 */
template <typename T>
void compare_errors(const T* leading_error, const T* absolute_errors_h, const std::size_t count_nodes) {
    if (leading_error == nullptr || absolute_errors_h == nullptr) throw std::invalid_argument("Input arrays cannot be null.");
    if (count_nodes <= 0) throw std::invalid_argument("Invalid count_nodes value.");
    // Вычисляем нормы для главного члена погрешности
    T* norms_leading_error = calculate_norms(leading_error, count_nodes);
    std::cout << std::string(61, '-') << "\n";
    std::cout << std::setw(18) << "Leading Error"
              << std::scientific << std::setprecision(6)
              << std::setw(15) << norms_leading_error[0]
              << std::setw(15) << norms_leading_error[1]
              << std::setw(15) << norms_leading_error[2] << "\n";

    std::cout << std::setw(18) << "Absolute Errors"
              << std::setw(15) << absolute_errors_h[0]
              << std::setw(15) << absolute_errors_h[1]
              << std::setw(15) << absolute_errors_h[2] << "\n";

    std::cout << std::string(61, '-') << "\n";

    delete[] norms_leading_error; // Освобождаем память
}
