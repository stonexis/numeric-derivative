#pragma once
#include <string>

namespace Task_const {
    /// Редактируемые параметры
    inline constexpr long double A = -4.0; ///Концы отрезка
    inline constexpr long double B = 4.0; ///Концы отрезка
    inline constexpr std::size_t M = 30; ///Количество узлов сетки
    inline constexpr std::size_t M_viz_ratio = 20; ///Количество раз, во сколько измельчить сетку, для отображения графика

    /// Нередактируемые параметры 
    inline constexpr std::size_t M_2 = M*2 - 1; /// Количество узлов сетки h/2
    inline const long double H = std::abs(B-A)/(M-1); ///Шаг равномерной сетки
    inline const long double STEP_H_2 = Task_const::H / 2; /// Шаг сетки h/2 
}

template <typename T>
std::tuple<T*, T*, T*> gen_grid_func_and_analytic_derivative(
                                    std::size_t& count_nodes_out,
                                    std::tuple<T*, T*, T*> x_func_derivative_rare=std::make_tuple(nullptr, nullptr, nullptr), 
                                    const std::size_t count_nodes_init=Task_const::M,
                                    const std::size_t ratio=1, 
                                    const T a=Task_const::A, const T b=Task_const::B
                                    );
template <typename T> 
T* gen_uniform_grid(
                const T step=Task_const::H,
                const std::size_t count_nodes=Task_const::M, 
                const T a=Task_const::A, 
                const T b=Task_const::B
                );
template <typename T>
std::pair<T*, T*> gen_runge_romberg(
                        const T* grid_derivative_more_freq, 
                        const T* grid_derivative_less_freq, 
                        const std::size_t count_nodes_great_freq=Task_const::M_2,
                        const std::size_t count_nodes_small_freq=Task_const::M,
                        const T step_greater_freq=Task_const::STEP_H_2,
                        const T step_smaller_freq=Task_const::H
                        );

template <typename T>
void write_to_file_arr(std::ofstream &out, const T *array, std::size_t length, std::string name_array);
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
                    );
template <typename T>
std::pair<T*,T*> calculate_norms(const T* analytical, const T* numerical, const std::size_t count_nodes);
template <typename T>
T* calculate_norms(const T* numerical, const std::size_t count_nodes);
template <typename T>
void print_error_table(
                const std::pair<T*,T*>& errors_h, 
                const std::pair<T*,T*>& errors_h_2, 
                const std::pair<T*,T*>& errors_runge, 
                const T* norms_leading_error
                );
#include "numerical_differentiation.tpp"