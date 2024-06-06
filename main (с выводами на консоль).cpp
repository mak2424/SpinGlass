/**
  Спиновое стекло.
  Модель Изинга.
  Открытые граничные условия.
**/

#include <iostream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <map>
#include <vector>

#define n 6
#define N n *n

#define cluster_size 3

#define MC_steps N * 10 // задаем количество Монте-Карло шагов

#define logs false

void J_glass_generator(int *J_h, int *J_v, int J_ferr);
void Ising_array_init(int *S);
void print_system(int *S);
void print_system(int *S, int *J_h, int *J_v);
void print_array(int *S, int size);
void print_system_arrow(int *S, int *J_h, int *J_v);
void print_cluster_arrow(int *cluster_array, int *S);
void print_cluster_arrow(int *cluster_state);
int E_system(int *S, int *J_h, int *J_v);
int E_spin(int spin, int *S, int *J_h, int *J_v);

// функция перебора всех состояний кластера
void cluster_bruteforce(int **cluster_states);
// функция определения всех кластеров
void get_square_clusters(int **clusters_array); // clusters_array[номер кластера][спины кластера]
// функция определения всех кластеров, их энергии и спинового избытка
void get_square_clusters(int **clusters_array, int *S, int *J_h, int *J_v);
// функция расчета энергии кластера
int E_cluster(int *cluster_array, int *S, int *J_h, int *J_v);
int E_cluster2(int *cluster_array, int *cluster_values, int *S, int *J_h, int *J_v);
// функция расчета спинового избытка кластера
int M_cluster(int *cluster_array, int *S);
int M_cluster2(int *cluster_array);

void HM(int *S, int *J_h, int *J_v, int **clusters_array, int num_of_clusters);

int main()
{
    int *J_h_host = new int[n * (n - 1)];
    int *J_v_host = new int[n * (n - 1)];
    int J_ferr = n * (n - 1); /// n * (n - 1);

    if (J_ferr > 2 * (n * (n - 1)))
    {
        std::cout << "to match J_ferr " << "\n";
        return 0;
    }

    J_glass_generator(J_h_host, J_v_host, J_ferr);

    unsigned int n_J = n * (n - 1);

    int *S = new int[N];
    Ising_array_init(S);
    // print_array(S);
    print_system_arrow(S, J_h_host, J_v_host);

    std::cout << "E_system = " << E_system(S, J_h_host, J_v_host) << std::endl;
    // std::cout << "E_spin[0] = " << E_spin(0, S, J_h_host, J_v_host) << std::endl;

    if (logs)
    {
        std::cout << "Coupling count: " << 2 * (n * (n - 1)) << std::endl;
        std::cout << std::endl;

        std::cout << "J_h_host: \n";
        for (unsigned int i = 0; i < n_J; i++)
            std::cout << J_h_host[i] << " ";

        std::cout << "\n\nJ_v_host: \n";
        for (unsigned int i = 0; i < n_J; i++)
            std::cout << J_v_host[i] << " ";
        std::cout << std::endl
                  << std::endl;
    }

    int n_exl = 0;
    for (int i = 0; i < cluster_size - 1; i++)
    {
        n_exl += (n - i) + (n - 1 - i);
    }
    int count_of_clusters = N - n_exl;
    std::cout << std::endl
              << "Count of clusters in the system = " << count_of_clusters << std::endl
              << std::endl;

    int count_of_spins_in_cluster = cluster_size * cluster_size;

    int **clusters_array = new int *[count_of_clusters];
    for (int count = 0; count < count_of_clusters; count++)
        clusters_array[count] = new int[count_of_spins_in_cluster];

    get_square_clusters(clusters_array);
    /// get_square_clusters(clusters_array, S, J_h_host, J_v_host);
    // std::cout << "E_cluster = " << E_cluster(clusters_array[0], S, J_h_host, J_v_host)<<std::endl;

    HM(S, J_h_host, J_v_host, clusters_array, count_of_clusters);

    delete[] S;
    delete[] J_h_host;
    delete[] J_v_host;

    for (int count = 0; count < count_of_clusters; count++)
        delete[] clusters_array[count];

    return 0;
} // main

// функция распределения обменных интегралов
void J_glass_generator(int *J_h, int *J_v, int J_ferr)
{
    for (auto i = 0; i < n * (n - 1); ++i)
    {
        J_h[i] = -1;
        J_v[i] = -1;
    }
    unsigned int n_J = n * (n - 1);
    std::random_device dev;
    std::mt19937 rand_gen(dev());
    rand_gen.seed(time(NULL));
    int *ferr_idx_h = new int[n_J];
    int *ferr_idx_v = new int[n_J];
    for (unsigned int i = 0; i < n_J; i++)
    {
        ferr_idx_h[i] = i;
        ferr_idx_v[i] = i;
    }
    unsigned int rem_h = n_J;
    unsigned int rem_v = n_J;
    for (auto i = 0; i < J_ferr; i++)
    {
        std::uniform_int_distribution<std::mt19937::result_type> rand_n(0, 2 * n_J - 1 - i);
        unsigned int idx = rand_n(rand_gen);
        if (idx < rem_h)
        {
            J_h[ferr_idx_h[idx]] = 1;
            rem_h--;
            for (auto j = idx; j < rem_h; j++)
            {
                ferr_idx_h[j] = ferr_idx_h[j + 1];
            }
        }
        else
        {
            J_v[ferr_idx_v[idx - rem_h]] = 1;
            rem_v--;
            for (auto j = idx - rem_h; j < rem_v; j++)
            {
                ferr_idx_v[j] = ferr_idx_v[j + 1];
            }
        }
    }
    delete[] ferr_idx_h;
    delete[] ferr_idx_v;
}

// инициализация системы спинов Изинга
void Ising_array_init(int *S)
{
    std::random_device dev;
    std::mt19937 rand_gen(dev());
    rand_gen.seed(time(NULL));
    std::uniform_int_distribution<std::mt19937::result_type> rand_n(0, 1);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (rand_n(rand_gen))
                S[i * n + j] = 1;
            else
                S[i * n + j] = -1; //-1;
        }
    }
}

// вывод спинов Изинга в виде +-1 на консоль
void print_system(int *S)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << S[i * n + j] << "\t";
        }
        std::cout << std::endl;
    }
}

// вывод спинов Изинга в виде +-1 на консоль
void print_array(int *S, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << S[i * size + j] << "\t";
        }
        std::cout << std::endl;
    }
}

// вывод спинов Изинга и обменных интегралов в виде +-1 на консоль
void print_system(int *S, int *J_h, int *J_v)
{
    int J_h_count = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << S[i * n + j];
            if (j < n - 1)
            {
                // std::cout << "\t" << J_h_count;
                // std::cout << "\t" << J_h[J_h_count];
                if (J_h[J_h_count] < 0)
                    std::cout << "\t" << "-" << "\t";
                else
                    std::cout << "\t" << "+" << "\t";
                J_h_count++;
            }
        }
        std::cout << std::endl
                  << std::endl;
        if (i < n - 1)
        {
            for (int j = 0; j < n; j++)
            {
                if (J_v[i * n + j] < 0)
                    std::cout << "-" << "\t\t";
                else
                    std::cout << "+" << "\t\t";
            }
        }
        std::cout << std::endl
                  << std::endl;
    }
}

// вывод спинов Изинга в виде стрелок и обменных интегралов на консоль
void print_system_arrow(int *S, int *J_h, int *J_v)
{
    int J_h_count = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (S[i * n + j] < 0)
                std::cout << (char)25;
            else
                std::cout << (char)24;
            if (j < n - 1)
            {
                // std::cout << "\t" << J_h_count;
                // std::cout << "\t" << J_h[J_h_count];
                if (J_h[J_h_count] < 0)
                    std::cout << " " << "-" << " ";
                else
                    std::cout << " " << "+" << " ";
                J_h_count++;
            }
        }

        std::cout << "\t\t";

        for (int j = 0; j < n; j++)
            std::cout << i * n + j << "\t";

        std::cout << std::endl;
        if (i < n - 1)
        {
            for (int j = 0; j < n; j++)
            {
                if (J_v[i * n + j] < 0)
                    std::cout << "-" << "   ";
                else
                    std::cout << "+" << "   ";
            }
        }
        std::cout << std::endl;
    }
}

// вывод спинов Изинга в виде стрелок на консоль
void print_cluster_arrow(int *cluster_state)
{
    for (int i = 0; i < cluster_size; i++)
    {
        for (int j = 0; j < cluster_size; j++)
        {
            if (cluster_state[i * cluster_size + j] < 0)
                std::cout << (char)25;
            else
                std::cout << (char)24;
            std::cout << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// вывод спинов Изинга в виде стрелок на консоль
void print_cluster_arrow(int *cluster_array, int *S)
{
    for (int i = 0; i < cluster_size; i++)
    {
        for (int j = 0; j < cluster_size; j++)
        {
            if (S[cluster_array[i * cluster_size + j]] < 0)
                std::cout << (char)25;
            else
                std::cout << (char)24;
            std::cout << "\t";
        }
        std::cout << "\t";
        for (int j = 0; j < cluster_size; j++)
        {
            std::cout << cluster_array[i * cluster_size + j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// расчет энергии системы с учетом открытых граничных условий
int E_system(int *S, int *J_h, int *J_v)
{
    if (logs)
        std::cout << "i\tj\t[i*n+j]\t[i*(n-1)+j]\t[i*n+j+1]\t[i*(n-1)+j+n]" << std::endl;
    int E_h = 0, E_v = 0, E;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < (n - 1); j++)
        {
            E_h += -J_h[i * (n - 1) + j] * S[i * n + j] * S[i * n + j + 1];
            E_v += -J_v[i * (n - 1) + j] * S[i * (n - 1) + j] * S[i * (n - 1) + j + n];
            if (logs)
                std::cout << i << "\t" << j << "\t" << i * n + j << "\t" << i * (n - 1) + j << "\t\t" << i * n + j + 1 << "\t\t" << i * (n - 1) + j + n << std::endl;
        }
    }
    if (logs)
        std::cout << std::endl;

    E = E_h + E_v;
    return E;
}

// расчет энергии одного спина с учетом открытых граничных условий
int E_spin(int spin, int *S, int *J_h, int *J_v)
{
    if (spin >= N)
    {
        std::cout << "Error: spin >= N" << std::endl;
        std::cout << " in file '" << __FILE__ << "' at line " << __LINE__ << std::endl;
        exit(EXIT_FAILURE);
    }

    int E_h = 0, E_v = 0, E;

    if ((spin + 1) % n)
        E_h += -J_h[spin - (spin / n)] * S[spin] * S[spin + 1]; // справа
    if (spin % n)
        E_h += -J_h[spin - 1 - ((spin - 1) / n)] * S[spin] * S[spin - 1]; // слева
    if ((spin / n + 1) % n)
        E_v += -J_v[spin] * S[spin] * S[spin + n]; // снизу
    if (spin / n)
        E_v += -J_v[spin - n] * S[spin] * S[spin - n]; // сверху

    if (logs)
        std::cout << std::endl;

    //    for(int i=0; i<N; i++){
    //        if((i+1)%n) std::cout << i << "\t" << i-(i/n) << std::endl;
    //    }

    E = E_h + E_v;
    return E;
}

// функция перебора всех состояний кластера
void cluster_bruteforce(int **cluster_states)
{
    std::cout << "Cluster size: " << cluster_size << "x" << cluster_size << std::endl
              << std::endl;

    int num_of_spins_in_cluster = cluster_size * cluster_size;       // количество спинов в кластере
    unsigned long long states_amount = 1 << num_of_spins_in_cluster; // количество конфигураций кластера

    // перебор конфигураций кластера побитовым сдвигом
    for (unsigned long long state_num = 0; state_num < states_amount; ++state_num)
    {
        /// std::cout << "State: " << state_num << " ( ";
        for (int spin_num_in_cluster = (num_of_spins_in_cluster - 1); spin_num_in_cluster >= 0; --spin_num_in_cluster)
        {
            if (state_num & (1 << spin_num_in_cluster))
                cluster_states[state_num][num_of_spins_in_cluster - 1 - spin_num_in_cluster] = 1;
            else
                cluster_states[state_num][num_of_spins_in_cluster - 1 - spin_num_in_cluster] = -1;

            /// std::cout << cluster_array[state_num][num_of_spins_in_cluster-1-spin_num_in_cluster] << " "; //вывод конфигурации на экран
        }
        /// std::cout<<")"<<std::endl;
    }
}

// функция определения всех кластеров
// clusters_array[номер кластера][спины кластера]
void get_square_clusters(int **clusters_array)
{
    // std::cout<<"Cluster size: "<<cluster_size<<"x"<<cluster_size<<std::endl<<std::endl;
    // std::cout<<"n_exl = "<<n_exl<<std::endl;

    int n_s;   // номер спина, от которого строится кластер
    int n_s_c; // номер спина, который попадает в кластер
    int n_c;   // порядковый номер кластера

    for (int i = 0; i < n - (cluster_size - 1); i++)
    {
        for (int j = 0; j < n - (cluster_size - 1); j++)
        {
            n_s = i * n + j;
            n_c = i * (n - (cluster_size - 1)) + j;

            for (int k = 0; k < cluster_size; k++)
            {
                for (int l = 0; l < cluster_size; l++)
                {
                    n_s_c = n_s + (k * n + l);
                    clusters_array[n_c][k * cluster_size + l] = n_s_c;
                    /// std::cout << clusters_array[n_c][k * cluster_size + l] << " ";
                }
                /// std::cout << std::endl;
            }
            /// std::cout << std::endl;
        }
    }
}

// функция определения всех кластеров, их энергии и спинового избытка
// clusters_array[номер кластера][спины кластера]
void get_square_clusters(int **clusters_array, int *S, int *J_h, int *J_v)
{
    // std::cout<<"Cluster size: "<<cluster_size<<"x"<<cluster_size<<std::endl<<std::endl;
    // std::cout<<"n_exl = "<<n_exl<<std::endl;

    int n_s;   // номер спина, от которого строится кластер
    int n_s_c; // номер спина, который попадает в кластер
    int n_c;   // порядковый номер кластера

    for (int i = 0; i < n - (cluster_size - 1); i++)
    {
        for (int j = 0; j < n - (cluster_size - 1); j++)
        {
            n_s = i * n + j;
            n_c = i * (n - (cluster_size - 1)) + j;
            std::cout << n_c << ":" << std::endl;

            for (int k = 0; k < cluster_size; k++)
            {
                for (int l = 0; l < cluster_size; l++)
                {
                    n_s_c = n_s + (k * n + l);
                    clusters_array[n_c][k * cluster_size + l] = n_s_c;
                    std::cout << clusters_array[n_c][k * cluster_size + l] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << "E_cluster = " << E_cluster(clusters_array[n_c], S, J_h, J_v) << std::endl;
            std::cout << "M_cluster = " << M_cluster(clusters_array[n_c], S) << std::endl;

            std::cout << std::endl;
        }
    }
}

// функция расчета энергии кластера
// cluster_array[спины кластера]
int E_cluster(int *cluster_array, int *S, int *J_h, int *J_v)
{
    // std::cout << "i\tj\t[i*n+j]\t[i*(n-1)+j]\t[i*n+j+1]\t[i*(n-1)+j+n]" << std::endl;
    int E_h = 0, E_v = 0, E, spin, cl_spin;
    for (int i = 0; i < cluster_size; i++)
    {
        for (int j = 0; j < (cluster_size); j++)
        {
            cl_spin = i * cluster_size + j; // порядковый номер спина в кластере
            spin = cluster_array[cl_spin];  // номер спина в системе
            if (i == 0)                     // верхняя грань
            {
                if (spin / n)
                    E_v += -J_v[spin - n] * S[spin] * S[spin - n]; // сверху
                if ((spin / n + 1) % n)
                    E_v += -J_v[spin] * S[spin] * S[spin + n]; // снизу
            }
            else
            {
                if ((spin / n + 1) % n)
                    E_v += -J_v[spin] * S[spin] * S[spin + n]; // снизу
            }

            if (j == 0) // левая грань
            {
                if (spin % n)
                    E_h += -J_h[spin - 1 - ((spin - 1) / n)] * S[spin] * S[spin - 1]; // слева
                if ((spin + 1) % n)
                    E_h += -J_h[spin - (spin / n)] * S[spin] * S[spin + 1]; // справа
            }
            else
            {
                if ((spin + 1) % n)
                    E_h += -J_h[spin - (spin / n)] * S[spin] * S[spin + 1]; // справа
            }

            /// E_h += -J_h[i*(cluster_size-1)+j] * S[cluster_array[i*cluster_size+j]] * S[cluster_array[i*cluster_size+j+1]];
            /// E_v += -J_v[i*(cluster_size-1)+j] * S[cluster_array[i*(cluster_size-1)+j]] * S[cluster_array[i*(cluster_size-1)+j+cluster_size]];
            // std::cout << i << "\t" << j << "\t" << i*cluster_size+j << "\t" << i*(cluster_size-1)+j << "\t\t" << i*cluster_size+j+1 << "\t\t" << i*(cluster_size-1)+j+cluster_size << std::endl;
        }
    }
    // std::cout<< std::endl;

    E = E_h + E_v;
    return E;
}

// функция расчета энергии кластера
// cluster_array[спины кластера]
int E_cluster2(int *cluster_array, int *cluster_state, int *S, int *J_h, int *J_v)
{
    // std::cout << "i\tj\t[i*n+j]\t[i*(n-1)+j]\t[i*n+j+1]\t[i*(n-1)+j+n]" << std::endl;
    int E_h = 0, E_v = 0, E, spin, cl_spin;
    int count_of_spins_in_cluster = cluster_size * cluster_size;
    for (int i = 0; i < cluster_size; i++)
    {
        for (int j = 0; j < cluster_size; j++)
        {
            cl_spin = i * cluster_size + j; // порядковый номер спина в кластере
            spin = cluster_array[cl_spin];  // порядковый номер спина в системе
            //std::cout << "[" << spin << "] " << cluster_state[cl_spin] << ", ";
            //system("pause");
            if (i == 0) // верхняя грань
            {
                if (spin / n)
                    E_v += -J_v[spin - n] * cluster_state[cl_spin] * S[spin - n]; // сверху
                if ((spin / n + 1) % n)
                    if (i == cluster_size - 1)
                        E_v += -J_v[spin] * cluster_state[cl_spin] * S[spin + n]; // снизу
                    else
                        E_v += -J_v[spin] * cluster_state[cl_spin] * cluster_state[cl_spin + cluster_size]; // снизу
            }
            else
            {
                if ((spin / n + 1) % n)
                    if (i == cluster_size - 1)
                        E_v += -J_v[spin] * cluster_state[cl_spin] * S[spin + n]; // снизу
                    else
                        E_v += -J_v[spin] * cluster_state[cl_spin] * cluster_state[cl_spin + cluster_size]; // снизу
            }

            if (j == 0) // левая грань
            {
                if (spin % n)
                    E_h += -J_h[spin - 1 - ((spin - 1) / n)] * cluster_state[cl_spin] * S[spin - 1]; // слева
                if ((spin + 1) % n)
                    if (j == cluster_size - 1)
                        E_h += -J_h[spin - (spin / n)] * cluster_state[cl_spin] * S[spin + 1]; // справа
                    else
                        E_h += -J_h[spin - (spin / n)] * cluster_state[cl_spin] * cluster_state[cl_spin + 1]; // справа
            }
            else
            {
                if ((spin + 1) % n)
                    if (j == cluster_size - 1)
                        E_h += -J_h[spin - (spin / n)] * cluster_state[cl_spin] * S[spin + 1]; // справа
                    else
                        E_h += -J_h[spin - (spin / n)] * cluster_state[cl_spin] * cluster_state[cl_spin + 1]; // справа
            }

            // E_h += -J_h[i*(cluster_size-1)+j] * S[cluster_array[i*cluster_size+j]] * S[cluster_array[i*cluster_size+j+1]];
            // E_v += -J_v[i*(cluster_size-1)+j] * S[cluster_array[i*(cluster_size-1)+j]] * S[cluster_array[i*(cluster_size-1)+j+cluster_size]];
            // std::cout << i << "\t" << j << "\t" << i*cluster_size+j << "\t" << i*(cluster_size-1)+j << "\t\t" << i*cluster_size+j+1 << "\t\t" << i*(cluster_size-1)+j+cluster_size << std::endl;
        }
    }

    E = E_h + E_v;

    //std::cout << "E_cl = " << E << std::endl;
    //std::cout << std::endl;
    return E;
}

// функция расчета спинового избытка кластера
int M_cluster(int *cluster_array, int *S)
{
    int M = 0, spin;
    for (int i = 0; i < cluster_size; i++)
    {
        for (int j = 0; j < (cluster_size); j++)
        {
            spin = cluster_array[i * cluster_size + j];
            M += S[spin];
            // std::cout << S[spin] << "\t" << std::endl;
        }
    }
    // std::cout << "M = " << M << std::endl;
    return M;
}

int M_cluster2(int *cluster_array)
{
    int M = 0;
    for (int i = 0; i < cluster_size; i++)
    {
        for (int j = 0; j < (cluster_size); j++)
        {
            M += cluster_array[i * cluster_size + j];
        }
    }
    return M;
}

/// void HM()
/// 1. Выбираем случайный кластер
/// 2. Вычисляем E и М кластера
/// 3. Делаем перебор всех состояний. Вычисляем E, M и p(E,M) каждой конфигурации.
///    p(E,M) = exp(-(E/T))/Z;
///    Z = SUM(exp(-(E/T)));
/// либо набора конфигураций с одинаковыми параметрами
///    p(E,M) = g*exp(-(E/T))/Z;
///    Z = SUM(g*exp(-(E/T)));
/// 4. Вычисляем термодинамику с учетом всех вариантов кластера
/// 5. Генерируем случайное число r от 0 до 1
/// 6. Если r<p(E,M) = g*exp(-(E/T))/SUM(g*exp(-(E/T))); то принимаем новую конфигурацию.
///    Иначе оставляем старую.
///    Для 2^𝑁𝑐 конфигураций поочередно суммируем их вероятности до тех пор,
///    пока выбранное на предыдущем шаге случайное число не попадет в полученный интервал (рисунок 2.1).
///    По факту суммируем вероятности по порядку, пока сумма не превысит r.
/// 7. Выбираем случайным образом конфигурацию из набора с одинаковыми параметрами.
///    Последней добавленной вероятности соответствует набор конфигураций,
///    из которого выбираем случайную конфигурацию кластера.
/// 8. Обновляем конфигурацию кластера.
/// 9. Обновляем средние термодинамические значения.

// гибридный метод алгоритм Метрополиса + полный перебор
// TODO: не готов
void HM(int *S, int *J_h, int *J_v, int **clusters_array, int num_of_clusters)
{
    std::cout << "\n***HM method***\n";
    srand(time(0));
    int cluster_number; // номер кластера
    int E_cl;           // энергия кластера
    int M_cl;           // спиновый избыток кластера

    int num_of_spins_in_cluster = cluster_size * cluster_size;       // количество спинов в кластере
    unsigned long long states_amount = 1 << num_of_spins_in_cluster; // количество конфигураций кластера

    int **cluster_states = new int *[states_amount];
    for (unsigned long long count = 0; count < states_amount; count++)
        cluster_states[count] = new int[num_of_spins_in_cluster];

    cluster_bruteforce(cluster_states); // Делаем перебор всех состояний кластера

    // TODO: нужна структура, чтобы объединить
    struct states_p
    {
        double p;
        std::vector<unsigned int> states;
    };

    // std::map<int, int> E_count;
    // std::map<int, std::map<int, unsigned int>> map_gEM;
    // std::map<int, std::map<int, double>> map_pEM;
    // std::map<int, std::map<int, std::vector<unsigned int>>> map_vecEM;

    std::map<int, std::map<int, states_p>> map_structEM;

    double r; // 0-1
    double f; // вероятность выбора конфигурации
    int rand_state;
    int E = E_system(S, J_h, J_v);

    std::cout << "E_system_start = " << E << std::endl;

    for (double T = 0.021; T > 0; T -= 0.01)
    {
        for (unsigned long long step = 1; step <= MC_steps; step++)
        {
            std::cout << "\nT = " << T << ", ";
            std::cout << "MC Step = " << step << "\n\n";
            // 1. Выбираем случайный кластер
            cluster_number = rand() % num_of_clusters;
            // cluster_number = 0;

            // 2. Вычисляем E и М кластера
            E_cl = E_cluster(clusters_array[cluster_number], S, J_h, J_v);
            M_cl = M_cluster(clusters_array[cluster_number], S);

            std::cout << "Cluster number: " << cluster_number << std::endl;
            std::cout << "Cluster E: " << E_cl << std::endl;
            std::cout << "Cluster M: " << M_cl << std::endl;
            print_cluster_arrow(clusters_array[cluster_number], S);
            // print_array(clusters_array[cluster_number], cluster_size);
            std::cout << std::endl;

            // 3. Вычисляем E, M и p(E,M) каждой конфигурации кластера.
            //    p(E,M) = exp(-(E/T))/Z;
            //    Z = SUM(exp(-(E/T)));
            // либо набора конфигураций с одинаковыми параметрами
            //    p(E,M) = g*exp(-(E/T))/Z;
            //    Z = SUM(g*exp(-(E/T)));

            double Z = 0; // статсумма
            double p;     // вероятность
            int flag = 0;
            int Emin = E_cl;

            // собираем g(E,M)
            for (unsigned long long cl_state = 0; cl_state < states_amount; cl_state++)
            {
                E_cl = E_cluster2(clusters_array[cluster_number], cluster_states[cl_state], S, J_h, J_v);
                M_cl = M_cluster2(cluster_states[cl_state]);
                // map_gEM[E_cl][M_cl]++;
                // map_vecEM[E_cl][M_cl].push_back(cl_state);
                map_structEM[E_cl][M_cl].states.push_back(cl_state); // TODO: изменить на сохранение энергии всей системы
                // std::cout << "E_cl: " << E_cl << std::endl;
                // std::cout << "M_cl: " << M_cl << std::endl;
                // std::cout << std::endl;
                if (E_cl < Emin)
                    Emin = E_cl;
            }

            // собираем статсумму
            for (auto it : map_structEM)
            {
                for (auto it2 : it.second)
                {
                    // std::cout << it.first << ", ";            // E
                    // std::cout << it2.first << ", ";           // M
                    // std::cout << it2.second.states.size() << std::endl; // g
                    p = it2.second.states.size() * exp(-((it.first-Emin) / T)); // p(E,M) = g*exp(-(E/T))/Z;
                    Z += p;                                              // Z = SUM(g*exp(-(E/T)));
                    map_structEM[it.first][it2.first].p = p;
                }
            }

            // double k = 0;
            //  вычисляем вероятность p(E,M)
            for (auto it : map_structEM)
            {
                for (auto it2 : it.second)
                {
                    // std::cout << it.first << ", ";         // E
                    // std::cout << it2.first << ", ";        // M
                    std::cout << "p(" << it.first << ", " << it2.first << ", g" << it2.second.states.size() << ") = " << it2.second.p / Z << ", "; // << std::endl;  // p(E,M)
                    map_structEM[it.first][it2.first].p = it2.second.p / Z;
                    // k += it2.second.p / Z;
                }
            }
            // std::cout << "SUM(p(E,M)) = " << k << std::endl;
            std::cout << std::endl;

            // TODO: 4. Вычисляем термодинамику с учетом всех вариантов кластера

            /////////////////////////////

            // 5. Генерируем случайное число r от 0 до 1
            r = (double)rand() / RAND_MAX;
            /// std::cout << "r = " << r << std::endl;

            // 6. Если r<p(E,M), то принимаем новую конфигурацию.
            //    Иначе оставляем старую.
            //    Для 2^𝑁𝑐 конфигураций поочередно суммируем их вероятности до тех пор,
            //    пока выбранное на предыдущем шаге случайное число не попадет в полученный интервал,
            //    (т.е. суммируем вероятности по порядку, пока сумма не превысит r).

            f = 0;
            flag = 0;

            for (auto it : map_structEM)
            {
                for (auto it2 : it.second)
                {
                    // std::cout << it.first << ", ";         // E
                    // std::cout << it2.first << ", ";        // M
                    // std::cout << it2.second.p << std::endl;  // p(E,M)

                    if (flag == 0)
                    {
                        if (r < (f + it2.second.p))
                        {
                            std::cout << std::endl
                                      << "p(" << it.first << ", " << it2.first << ") = " << it2.second.p << std::endl;
                            f += it2.second.p;
                            flag = 1;
                            std::cout << r << " < " << f << std::endl;
                            /// std::cout << "states: " << it2.second.states.size() << std::endl;

                            /// 7. Выбираем случайным образом конфигурацию из набора с одинаковыми параметрами.
                            ///    Последней добавленной вероятности соответствует набор конфигураций,
                            ///    из которого выбираем случайную конфигурацию кластера.
                            rand_state = rand() % it2.second.states.size();
                            /// std::cout << "rand_state: " << rand_state << std::endl;
                            rand_state = it2.second.states[rand_state];
                            /// std::cout << "rand_state: " << rand_state << std::endl;
                            /*
                            for (auto it_vec : it2.second.states)
                            {
                                std::cout << it_vec << ", ";
                            }
    */
                            /// std::cout << std::endl;

                            print_cluster_arrow(cluster_states[rand_state]);
                            /// std::cout << std::endl;

                            /// 8. Обновляем конфигурацию кластера.
                            for (int i = 0; i < num_of_spins_in_cluster; i++)
                            {
                                S[clusters_array[cluster_number][i]] = cluster_states[rand_state][i];
                                /// std::cout << cluster_states[rand_state][i] << ", ";
                            }

                            std::cout << "E_system = " << E_system(S, J_h, J_v) << std::endl;

                            print_system_arrow(S, J_h, J_v);

                            system("pause");
                        }
                        else
                        {
                            f += it2.second.p;
                        }
                    }
                }
            }

            /// print_array_arrow(S, J_h, J_v);

            /// 9. Обновляем средние термодинамические значения.
            /// 10. Ищем конфигурацию с минимальной энергией.

            // map_vecEM.erase(map_vecEM.begin(), map_vecEM.end());
            // map_pEM.erase(map_pEM.begin(), map_pEM.end());
            map_structEM.erase(map_structEM.begin(), map_structEM.end());
            //system("pause");
        } // MC

    } // T

    for (int count = 0; count < num_of_spins_in_cluster; count++)
        delete[] cluster_states[count];
}

/*
void HM()
{
    ofstream C_data("C.txt"); //файл теплоемкости
    ofstream hi_data("hi.txt"); //файл магнитной восприимчивости
    ofstream E_data("E.txt"); //файл энергии

    double E_aver = 0, E2_aver = 0; //средняя энергия и ее квадрат
    double Mpr_aver = 0, Mpr2_aver = 0; //средняя проекция намагниченности и ее квадрат

    unsigned long long set_of_states, total_set_of_states=0;
    bool flag=0;

    ///uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение

    int rand_state=0; //случайная конфигурация с одинаковыми параметрами EM

    double Z = 0; //статсумма
    double r; //число от 0 до 1
    double E_aver_i = 0;
    double E2_aver_i = 0;
    double M_aver_i = 0;
    double M2_aver_i = 0;

    //интервалы вероятностей
    double interval=0;

    unsigned long long spin_num; //случайный спин

    ///uniform_int_distribution<int> distribution_int(0,(N*N)-1);

    //проход блока по системе
    for(double T = Tmax; T>=Tmin; T-=Tstep)
    //for(float T = Tmin; T<Tmax; T+=Tstep) //цикл по температуре
    {
        T = round(T*10000)/10000.;
        ///if(T>2.1 && T<2.5) Tstep=0.01;
        ///else Tstep=0.1;

        //T=2.269;//------------------------------------------------------------------------убрать

        cout<<"T = "<<T<<endl<<"----------"<<endl;

        E_aver = 0;
        E2_aver = 0;
        Mpr_aver = 0;
        Mpr2_aver = 0;


        for(unsigned long long MCS = 1; MCS<=NumMC; ++MCS)
        {
            spin_num = rand()%(N*N);

            ///ПЕРЕДЕЛАЛ!!! ПРОВЕРИТЬ!!!
            /// unsigned long long array_to_dec(int *array_1D) //Оформить как функцию
            unsigned long long boundary_dec=0; //конфигурация границы в десятичной системе
            int bit;
            for (int ii=0; ii<num_of_spins_in_boundaries; ++ii)
            {
                if (array_1D[boundary_coordinates_for_spins[spin_num][ii]]==-1)
                    bit = 0;
                else
                    bit = 1;

                boundary_dec = (boundary_dec | bit);

                if (ii<num_of_spins_in_boundaries-1)
                    boundary_dec = boundary_dec<<1;
            }
            ///cout<<"boundary_dec = "<<boundary_dec<<endl;

            unsigned long long block_dec=0; //конфигурация блока в десятичной системе
            for (int ii=0; ii<num_of_spins_in_block; ++ii)
            {
                if (array_1D[block_coordinates_for_spins[spin_num][ii]]==-1)
                    bit = 0;
                else
                    bit = 1;

                block_dec = (block_dec | bit);

                if (ii<num_of_spins_in_block-1)
                    block_dec = block_dec<<1;
            }

            unsigned long long num_of_unique_EM_in_boundary = boundaries_EM_blocks[boundary_dec].size(); //g
            double P[num_of_unique_EM_in_boundary]; //массив вероятностей энергий

            ///cout<<"E0 = "<<boundaries_blocks_EM[boundary_dec][block_dec].E<<", M0 = "<<boundaries_blocks_EM[boundary_dec][block_dec].M<<endl;
            ///cout << "Emin = " << Emin_array[boundary_dec] << endl;

            old_E-= boundaries_blocks_EM[boundary_dec][block_dec].E;
            My-= boundaries_blocks_EM[boundary_dec][block_dec].M;

            Z = 0; //статсумма
            E_aver_i = 0;
            E2_aver_i = 0;
            M_aver_i = 0;
            M2_aver_i = 0;

            //считаем статсумму
            for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
            {
                //вырождение * exp
                P[set_of_states] = boundaries_EM_blocks[boundary_dec][set_of_states].state_array.size()*
                        exp((double)(-(boundaries_EM_blocks[boundary_dec][set_of_states].E-Emin_array[boundary_dec])/T));
                Z += P[set_of_states];
                ///cout<<exp((double)(-(boundaries_EM_blocks[boundary_dec][set_of_states].E-Emin)/T))<<endl;
                ///
                /// ПЕРЕНЕСТИ СЮДА РАСЧЕТ СРЕДНИХ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            }

            ///system("pause");

            ///cout << "Z = " << Z << endl;

            ///r = distribution_real(generator);
            r = ((double) rand() / (RAND_MAX));

            ///cout<<"r= "<<r<<endl;

            //перебираем интервалы вероятностей
            interval=0;
            flag=0;
            //double sum=0;

            //считаем вероятности
            for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
            {
                P[set_of_states] /= Z;
                //sum+=P[set_of_states];
                ///cout << "P"<<set_of_states<<"= " << P[set_of_states] << endl;
                E_aver_i += P[set_of_states]*(old_E+boundaries_EM_blocks[boundary_dec][set_of_states].E);
                M_aver_i += P[set_of_states]*abs(My+boundaries_EM_blocks[boundary_dec][set_of_states].M);

                E2_aver_i += P[set_of_states]*(old_E+boundaries_EM_blocks[boundary_dec][set_of_states].E)*(old_E+boundaries_EM_blocks[boundary_dec][set_of_states].E);
                M2_aver_i += P[set_of_states]*(My+boundaries_EM_blocks[boundary_dec][set_of_states].M)*(My+boundaries_EM_blocks[boundary_dec][set_of_states].M);

                interval+=P[set_of_states];
                if((r<interval || is_equal(r,interval)) && flag==0)
                {
                    total_set_of_states = set_of_states;
                    flag = 1;
                }
            }
            ///cout<<"E_aver_i = "<<E_aver_i<<endl;
            //cout<<"P_sum = " << sum << endl;

            My += boundaries_EM_blocks[boundary_dec][total_set_of_states].M;
            old_E += boundaries_EM_blocks[boundary_dec][total_set_of_states].E;


            //средние
            //E_aver += (old_E - E_aver) / (double)MCS;
            E_aver += (E_aver_i - E_aver) / (double)MCS;
            //E2_aver += (old_E*old_E - E2_aver) / (double)MCS;
            //E2_aver += (E_aver_i*E_aver_i - E2_aver) / (double)MCS;
            E2_aver += (E2_aver_i - E2_aver) / (double)MCS;

            Mpr_aver += (M_aver_i - Mpr_aver) / (double)MCS;
            //Mpr_aver += (abs(M_aver_i) - Mpr_aver) / (double)MCS;
            //Mpr2_aver += (M_aver_i*M_aver_i - Mpr2_aver) / (double)MCS;
            Mpr2_aver += (M2_aver_i - Mpr2_aver) / (double)MCS;

            ///cout<<"E_aver = "<<E_aver<<endl;
            ///system("pause");

            ///cout<<"r<="<<interval<<endl;

            ///My += boundaries_EM_blocks[boundary_dec][total_set_of_states].M;
            ///old_E += boundaries_EM_blocks[boundary_dec][total_set_of_states].E;

            //выбираем случайную конфигурацию с одинаковыми параметрами
            ///uniform_int_distribution<int> distribution_int2(0,boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array.size()-1);
            ///rand_state = distribution_int2(generator);
            if (boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array.size()==0){
                cout<<"boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array.size() = "<<boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array.size()<<endl;
                system("pause");
            }
            rand_state = rand() % boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array.size();

            //запоминаем конфигурацию
            rand_state = boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array[rand_state];

            //переводим в двоичную систему
            ///cout<<"rand_state = "<<rand_state<<endl;

            ///ПЕРЕДЕЛАЛ!!! ПРОВЕРИТЬ!!!
            count=0;
            for(bit=num_of_spins_in_block-1; bit>=0; --bit)
            {
                array_1D[block_coordinates_for_spins[spin_num][count]] = (1 & rand_state >> bit);

                if(array_1D[block_coordinates_for_spins[spin_num][count]]==0)
                    array_1D[block_coordinates_for_spins[spin_num][count]] = -1;

                ///cout << block_temp[count]<< " ";
                count++;
            }
            ///cout<<endl;

        }//здесь заканчивается проход блока по системе
        //break;//------------------------------------------------------------


        //cout << "Mpr_aver = " << Mpr_aver << endl << endl;

        cout << "E_aver = " << E_aver << endl;

        double C = ((E2_aver - E_aver*E_aver)/(T*T))/(N*N);
        cout << "C = " << C << endl;
        C_data << T << "\t" << C << endl;
        //if(C1>C) Tstep /= 1.1;
        //else Tstep *= 1.1;
        //C1=C;

        double hi = ((Mpr2_aver - Mpr_aver*Mpr_aver)/T) / (double)(N*N);
        cout << "hi = " << hi << endl << endl;
        hi_data << T << "\t" << hi << endl;

        E_data << T << "\t" << E_aver << endl;

        ///system("pause");


        if(old_E < (-2*N*N) || old_E > (2*N*N)) cout<<"ERROR2 E = "<<old_E<<endl;
    }//здесь заканчивается перебор температур

    C_data.close();
    hi_data.close();
    E_data.close();

    ///cout << "\n" << (double)clock()/CLOCKS_PER_SEC << " sec.\n";
    ///system("pause");
}
*/
