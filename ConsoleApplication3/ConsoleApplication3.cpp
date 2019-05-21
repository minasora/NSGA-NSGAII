// ConsoleApplication3.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 用NSGA解决DEB函数优化
// 采用浮点数编码,模拟二进制交叉,多项式变异，锦标赛选择
// Huazhong University of Science and Technology 

#include "pch.h"
#include <iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<vector>
#include<algorithm>
using namespace std;
/********************常量*******************/
const int Virtual_fitness = 1;//指定第一级虚拟适应度为1，随后每级依次递减0.1
const double alpha = 1; //指定的共享半径
const int n = 20;//种群个体个数
const int gen = 10000;//迭代次数
const double high = 1.0;
const double low = 0.0; //解空间的范围;
const double deta_m = 1.0;//分布指数，DEB推荐为1.0,可自己调整
const double PI = 3.1415926;
/********************变量*******************/
struct Solution//解
{
	double x1;
	double x2;//两个决策变量
	double y1;
	double y2;//两个寻优函数
	double fitness;//解的虚拟适应度            
	bool if_rank;//是否已经排序
};
vector<vector<Solution>> dominate_rank;//级别
vector<Solution> population;//解的集合


/********************函数*******************/
double random_number(int x, int y)//返回[x,y]间的随机数
{
	return x + (rand() / RAND_MAX)*(y - x);
}
int rand_0_1()//随机生成0和1
{
	return rand() % 2;
}
void calculate_f(Solution &a)//计算目标函数值
{
	a.y1 = a.x1;
	a.y2 = (1 + 10 * a.x2)*(1 - (a.x1 / 1 + 10 * a.x2)*(a.x1 / 1 + 10 * a.x2) - a.x1 / 1 + 10 * a.x2*sin(8 * PI * a.x1));//DEB函数
}
void init()//初始化
{
	for (int i = 1; i <= n; i++)
	{
		Solution k;
		k.x1 = random_number(0, 1);
		k.x2 = random_number(0, 1);
		calculate_f(k);
		population.push_back(k);
	}
}
bool is_dominate(Solution a, Solution b)//判断两个解的支配关系
{
	if (a.y1 >= b.y1 && a.y2 >= b.y2)return true;
	else return false;
}
double distance_count(Solution a, Solution b)//计算两个解的欧式距离
{
	return sqrt((a.x1 - b.x1)*(a.x1 - b.x1)*(a.x2 - b.x2)*(a.x2 - b.x2));
}
double calculate_sharing(Solution a, Solution b)//计算一个解相对另一个解的共享度
{
	double sh;
	if (distance_count(a, b) > alpha)return sh = 0;//距离大于共享半径直接为0;
	else return sh = (1 - distance_count(a, b) / alpha);
}
void niche_sharing(vector<Solution> &S,int rank)//共享小生境
{
	for (int i = 0; i <= (S.size() - 1); i++)
	{
		double ans_sh = 0;//同级别共享适应度的和
		for (int j = 0; j <= (S.size() - 1); j++)
		{
			ans_sh += calculate_sharing(S[i], S[j]);
		}
		S[i].fitness = S[i].fitness / ans_sh;//重新计算每个个体的适应度，相近的收敛至小生境
	}
}
void Non_dominate_sort()//非支配排序
{
	int tmp = 1;//等级个数
	int num = 1;//判断是否完全排序
	while (true)
	{
		for (int i = 1; i <= n; i++)
		{
			if (population[i].if_rank)continue;
			for (int j = 1; j <= n; j++)
			{
				if (population[i].if_rank)continue;//发现定好等级的直接跳过
				if (!is_dominate(population[i], population[j]))//假如有非支配解直接跳出循环
				{
					break;
				}
				population[i].if_rank = true;
				num++;
				dominate_rank[tmp].push_back(population[i]);
			}
		}
		if (num == n)break;//排好序的数量和种群数量个数相等就跳出
		tmp++;
	}
}
void tournament_select()//每次随机选择两个个体优者保留下一代
{
	int i1 = rand() % dominate_rank.size();
	int j1 = rand() % dominate_rank[i1].size();
	int i2 = rand() % dominate_rank.size();
	int j2 = rand() % dominate_rank[i2].size();
	if (dominate_rank[i1][j1].fitness > dominate_rank[i2][j2].fitness)
	{
		population.push_back(dominate_rank[i1][j1]);
	}
	else
	{
		population.push_back(dominate_rank[i2][j2]);
	}
}
void crossever(Solution &a, Solution &b)//模拟二进制交叉
{
	double u = random_number(0, 1);//根据[0,1]间随机值决定了方程
	double deta_j;
	if (u < 0.5) 
	{
		deta_j = pow((2 * u), 1 / (deta_m + 1));
	}
	else
	{
		deta_j = pow(1 / (2 * (1 - u)), 1 / (deta_m + 1));
	}
	double child1x1 = 0.5*(((1 + deta_j)*a.x1 + (1 - deta_j)*b.x1));
	double child2x1 = 0.5*(((1 - deta_j)*a.x1 + (1 + deta_j)*b.x1));
	double child1x2 = 0.5*(((1 + deta_j)*a.x2 + (1 - deta_j)*b.x2));
	double child2x2 = 0.5*(((1 - deta_j)*a.x2 + (1 + deta_j)*b.x2));
	a.x1 = child1x1;
	a.x2 = child1x2;
	b.x1 = child2x1;
	b.x2 = child2x2;
}
void mutation(Solution &a)//多项式突变
{
	double u = random_number(0, 1);
	double delta1 = (a.x1 - low) / (high - low);
	double delta2 = (high - a.x1) / (high - low);
	double lambaq;
	if (u <= 0.5)
	{
		lambaq = pow(pow((2 * u + (1 - 2u)*(1 - delta1)), 1 + deta_m),1/(deta_m+1))-1.0;
	}
	else
	{
		lambaq = 1.0 - pow(pow((2 * u + (1 - 2u)*(1 - delta1)), 1 + deta_m), 1 / (deta_m + 1));
	}
	a.x1 += lambaq * (high - low);
    delta1 = (a.x2 - low) / (high - low);
	delta2 = (high - a.x2) / (high - low);
    lambaq;
	if (u <= 0.5)
	{
		lambaq = pow(pow((2 * u + (1 - 2u)*(1 - delta1)), 1 + deta_m), 1 / (deta_m + 1)) - 1.0;
	}
	else
	{
		lambaq = 1.0 - pow(pow((2 * u + (1 - 2u)*(1 - delta1)), 1 + deta_m), 1 / (deta_m + 1));
	}
	a.x2 += lambaq * (high - low);
}
int main()
{
	int gennumbers = 0;
	init();
	while (gennumbers<gen)
	{
		Non_dominate_sort();
		for (int i = 0; i <= (dominate_rank.size() - 1); i++)
		{
			for (int j = 0; j <= (dominate_rank.size() - 1); j++)
			{
				dominate_rank[i][j].fitness = Virtual_fitness - i * 0.1;//为每层指定不同的虚拟适应度;
			}
			niche_sharing(dominate_rank[i], i);//每层都共享适应度
		}
		population.clear();
		while (population.size() != n)
		{
			tournament_select();//选择
		}
		for (int i = 1; i <= n; i++)
		{
			int j = rand() % n;
			int k = rand() % n;
			crossever(population[j], population[k]);//随机交叉
		}
		for (int i = 1; i <= n; i++)
		{
			int j = rand() % n;
			mutation(population[j]);
		}
		for (int i = 1; i <= n; i++)
			calculate_f(population[i]);
		gennumbers++;
	}

}
