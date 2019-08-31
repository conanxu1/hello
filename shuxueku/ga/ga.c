/*
 * 遗传算法求解函数的最优值问题
 * 参考自《MATLAB智能算法30个案例分析》
 * 本例的寻优函数为:f = -5*sin(x1)*sin(x2)*sin(x3)*sin(x4)*sin(x5) - sin(5*x1)*sin(5*x2)*sin(5*x3)*sin(5*x4)*sin(5*x5) + 8
 * 其中x1,x2,x3,x4,x5是0~0.9*PI之间的实数。该函数的最小值为2,当x1,x2,x3,x4,x5都取PI/2时得到。
 * update in 16/12/3
 * author:Lyrichu
 * email:919987476@qq.com
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define PI 3.1415926 //圆周率
#define sizepop 50 // 种群数目
#define maxgen 500 // 进化代数
#define pcross 0.6 // 交叉概率
#define pmutation 0.1 // 变异概率
#define lenchrom 5 // 染色体长度
#define bound_down 0 // 变量下界，这里因为都相同，所以就用单个值去代替了，如果每个变量上下界不同，也许需要用数组定义
#define bound_up (0.9*3.1415926) // 上界
double chrom[sizepop][lenchrom]; // 种群数组
double fitness[sizepop]; //种群每个个体的适应度
double fitness_prob[sizepop]; // 每个个体被选中的概率(使用轮盘赌法确定)
double bestfitness[maxgen]; // 每一代最优值
double gbest_pos[lenchrom]; // 取最优值的位置
double average_best[maxgen+1]; // 每一代平均最优值
double gbest; // 所有进化中的最优值
int gbest_index; // 取得最优值的迭代次数序号


// 目标函数
double fit_func(double * arr)
{
    double x1 = *arr;
    double x2 = *(arr+1);
    double x3 = *(arr+2);
    double x4 = *(arr+3);
    double x5 = *(arr+4);
    double func_result = -5*sin(x1)*sin(x2)*sin(x3)*sin(x4)*sin(x5) - sin(5*x1)*sin(5*x2)*sin(5*x3)*sin(5*x4)*sin(5*x5) + 8;
    return func_result;
}

// 求和函数
double sum(double * fitness)
{
    double sum_fit = 0;
    for(int i=0;i<sizepop;i++)
        sum_fit += *(fitness+i);
    return sum_fit;
}
// 求最小值函数
double * min(double * fitness)
{
    double min_fit = *fitness;
    double min_index = 0;
    static double arr[2];
    for(int i=1;i<sizepop;i++)
    {
        if(*(fitness+i) < min_fit)
        {
            min_fit = *(fitness+i);
            min_index = i;
        }
    }
    arr[0] = min_index;
    arr[1] = min_fit;
    return arr;
}
// 种群初始化
void init_chrom()
{
    for(int i=0;i<sizepop;i++)
    {
        for(int j=0;j<lenchrom;j++)
        {
            chrom[i][j] = bound_up*(((double)rand())/RAND_MAX);
        }
        fitness[i] = fit_func(chrom[i]); // 初始化适应度
    }
}
// 选择操作
void Select(double chrom[sizepop][lenchrom])
{
    int index[sizepop];
    for(int i=0;i<sizepop;i++)
    {
        double * arr = chrom[i];
        fitness[i] = 1/(fit_func(arr)); // 本例是求最小值，适应度为目标函数的倒数，即函数值越小，适应度越大
    }
    double sum_fitness = 0;
    for(int i=0;i<sizepop;i++)
    {
        sum_fitness += fitness[i]; // 适应度求和
    }
    for(int i=0;i<sizepop;i++)
    {
        fitness_prob[i] = fitness[i]/sum_fitness;
    }
    for(int i=0;i<sizepop;i++)
    {
        fitness[i] = 1/fitness[i]; // 恢复函数值
    }
    for(int i=0;i<sizepop;i++) // sizepop 次轮盘赌
    {
        double pick = ((double)rand())/RAND_MAX;
        while(pick < 0.0001)
            pick = ((double)rand())/RAND_MAX;
        for(int j=0;j<sizepop;j++)
        {
            pick = pick - fitness_prob[j];
            if(pick<=0)
            {
                index[i] = j;
                break;
            }

        }
        
    }
    for(int i=0;i<sizepop;i++)
    {
        for(int j=0;j<lenchrom;j++)
        {
            chrom[i][j] = chrom[index[i]][j];
        }
        fitness[i] = fitness[index[i]];
    }
}

// 交叉操作
void Cross(double chrom[sizepop][lenchrom])
{
    for(int i=0;i<sizepop;i++)
    {
        // 随机选择两个染色体进行交叉
        double pick1 = ((double)rand())/RAND_MAX;
        double pick2 = ((double)rand())/RAND_MAX;
        int choice1 = (int)(pick1*sizepop);// 第一个随机选取的染色体序号
        int choice2 = (int)(pick2*sizepop);// 第二个染色体序号
        while(choice1 > sizepop-1)
        {
            pick1 = ((double)rand())/RAND_MAX;
            choice1 = (int)(pick1*sizepop); //防止选取位置过界
        }
        while(choice2 > sizepop-1)
        {
            pick2 = ((double)rand())/RAND_MAX;
            choice2 = (int)(pick2*sizepop);
        }
        double pick = ((double)rand())/RAND_MAX;// 用于判断是否进行交叉操作
        if(pick>pcross)
            continue;
        int flag = 0; // 判断交叉是否有效的标记
        while(flag == 0)
        {
            double pick = ((double)rand())/RAND_MAX;
            int pos = (int)(pick*lenchrom);
            while(pos > lenchrom-1)
            {
                double pick = ((double)rand())/RAND_MAX;
                int pos = (int)(pick*lenchrom); // 处理越界
            }
            // 交叉开始
            double r = ((double)rand())/RAND_MAX;
            double v1 = chrom[choice1][pos];
            double v2 = chrom[choice2][pos];
            chrom[choice1][pos] = r*v2 + (1-r)*v1;
            chrom[choice2][pos] = r*v1 + (1-r)*v2; // 交叉结束
            if(chrom[choice1][pos] >=bound_down && chrom[choice1][pos]<=bound_up && chrom[choice2][pos] >=bound_down && chrom[choice2][pos]<=bound_up)
                flag = 1; // 交叉有效，退出交叉，否则继续下一次交叉

        }

    }
}

// 变异操作
void Mutation(double chrom[sizepop][lenchrom])
{
    for(int i=0;i<sizepop;i++)
    {
        double pick = ((double)rand())/RAND_MAX; // 随机选择一个染色体进行变异
        int choice = (int)(pick*sizepop);
        while(choice > sizepop-1)
        {
            pick = ((double)rand())/RAND_MAX;
            int choice = (int)(pick*sizepop);  // 处理下标越界
        }
        pick = ((double)rand())/RAND_MAX; // 用于决定该轮是否进行变异
        if(pick>pmutation)
            continue;
        pick = ((double)rand())/RAND_MAX;
        int pos = (int)(pick*lenchrom);
        while(pos > lenchrom-1)
        {
            pick = ((double)rand())/RAND_MAX;
            pos = (int)(pick*lenchrom);
        }
        double v = chrom[i][pos];
        double v1 = v - bound_up;
        double v2 = bound_down - v;
        double r = ((double)rand())/RAND_MAX;
        double r1 = ((double)rand())/RAND_MAX;
        if(r >= 0.5)
            chrom[i][pos] = v - v1*r1*(1-((double)i)/maxgen)*(1-((double)i)/maxgen);
        else
            chrom[i][pos] = v + v2*r1*(1-((double)i)/maxgen)*(1-((double)i)/maxgen); // 注意这里写法是不会越界的，所以只用一次变异就可以了
    }
}

// 主函数
int main(void)
{
    time_t start,finish;
    start = clock(); // 程序开始计时
    srand((unsigned)time(NULL)); // 初始化随机数种子
    init_chrom(); // 初始化种群
    double * best_fit_index = min(fitness);
    int best_index =(int)(*best_fit_index);
    gbest = *(best_fit_index+1); // 最优值
    gbest_index = 0;
    average_best[0] = sum(fitness)/sizepop; //初始平均最优值
    for(int i=0;i<lenchrom;i++)
        gbest_pos[i] = chrom[best_index][i];
    // 进化开始
    for(int i=0;i<maxgen;i++)
    {
        Select(chrom); // 选择
        Cross(chrom); //交叉
        Mutation(chrom); //变异
        for(int j=0;j<sizepop;j++)
        {
            fitness[j] = fit_func(chrom[j]);
        }
        double sum_fit = sum(fitness);
        average_best[i+1] = sum_fit/sizepop; // 平均最优值
        double * arr = min(fitness);
        double new_best = *(arr+1);
        int new_best_index = (int)(*arr); // 新最优值序号
        if(new_best < gbest)
        {
            gbest = new_best;
            for(int j=0;j<lenchrom;j++)
            {
                gbest_pos[j] = chrom[new_best_index][j];
            }
            gbest_index = i+1;
        }

    }
    finish = clock(); // 程序计算结束
    double duration = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("程序计算耗时:%lf秒\n.",duration);
    printf("遗传算法进化了%d次，最优值为:%lf,最优值在第%d代取得,此代的平均最优值为%lf.\n",maxgen,gbest,gbest_index,average_best[gbest_index]);
    printf("取得最优值的地方为(%lf,%lf,%lf,%lf,%lf).\n",gbest_pos[0],gbest_pos[1],gbest_pos[2],gbest_pos[3],gbest_pos[4]);
    return 0;
}