import json
import seaborn as sns
import matplotlib.pyplot as plt


def plot_graph_seaborn(filename):
    # Чтение JSON данных из файла
    with open(filename, 'r') as file:
        data = json.load(file)

    # Преобразование двумерных массивов в одномерные
    grid_M_viz = data['grid_M_viz']
    grid_h = data['grid_h']
    grid_h_2 = data['grid_h_2']
    derivative_analytics = data['derivative_analytics']
    derivative_in_h = data['derivative_in_h']
    derivative_in_h_2 = data['derivative_in_h_2']
    updated_runge = data['updated_runge']

    # Установка темы и стиля
    sns.set_theme(style="whitegrid", palette="muted", font="serif", font_scale=1.2)

    # Построение графика
    plt.figure(figsize=(14, 8))

    # Построение линий для исходной функции и аппроксимации
    sns.lineplot(x=grid_M_viz, y=derivative_analytics, label='analyt', linewidth=2, zorder=1)
    sns.scatterplot(x=grid_h, y=derivative_in_h, s=50, label='h', marker='o', zorder=2, color='purple')
    sns.scatterplot(x=grid_h_2, y=derivative_in_h_2, s=40, label='h_2', marker='o', zorder=3, color='red')
    sns.scatterplot(x=grid_h_2, y=updated_runge, s=30, label='runge', marker='o', zorder=4, color='green')

    # Настройки графика
    plt.xlabel("X", fontsize=14)
    plt.ylabel("Y", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Настройка сетки
    plt.grid(which='major', linestyle='-', linewidth=0.5)  # Основная сетка
    plt.grid(which='minor', linestyle=':', linewidth=0.5)  # Мелкая сетка
    plt.minorticks_on()  # Включение дополнительной сетки
    plt.tight_layout()

    # Отображение графика
    plt.show()

# Вызов функции
plot_graph_seaborn("data.json")


