#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace std;

const int N = 10;				//粒子个数
const int DIM = 2;				//问题维数
const int MAX_STEP = 10000000;	//最大迭代步数
const double Ms = 0.2;			//Nose-Hover热库参数
const double Temperature =1.0;	//Nose-Hover热库温度
const bool isNoseHover = true;	//Nose-Hover热库开关

const double PI = 3.1415926;

double* pos;				//存储位置
double*	pre_pos;			//存储上一步位置
double* velocity;			//存储速度（暂未使用）
double* force;				//存储受力
double* potential;
double* kinetic_energy;
double max_force;			//存储最大受力
double eta = -1.0;
double time_step;	//时间步长
double time=0.0;	//统计当前物理时间
int step;

void initialize(){
	//此处设置初始位置
	//初始位置按照圆周均匀分布
	for (int i = 0; i < 10; i++) {
		pos[i*DIM + 0] = cos(2*PI * ((double)i) / N);
		pos[i*DIM + 1] = sin(2*PI * ((double)i) / N);
	}

	for(int i=0;i<DIM*N; i++)
		pre_pos[i] = pos[i];	

	for(int i=0;i<DIM*N; i++)
		*(velocity+i) = 0.0;
}

void compute_force(){
	/*
	  totalForce 两个粒子ij之间的相互作用力合力,吸引为负
	  dist 两个粒子之间的距离
	  */
	for(int i=0; i<N*DIM; i++){
		force[i] = 0.0;				//将每个粒子的受力清零
		max_force = 0.0;			//存储最大合力（用于变步长计算）
	}
	double  r, x, y, dist, dx, dy, totalForce;
	for(int i=0;i<N;i++){
		x = pos[i*DIM + 0];
		y = pos[i*DIM + 1];
		r = x*x + y*y;
		force[i*DIM + 0] += -x;		//外势力作用
		force[i*DIM + 1] += -y;		//外势力作用
		for(int j=i+1;j<N;j++){
			dx = pos[i*DIM + 0] - pos[j*DIM + 0];
			dy = pos[i*DIM + 1] - pos[j*DIM + 1];
			dist = sqrt(dx*dx+dy*dy);
			totalForce = (48*pow(dist,-13.0)-24*pow(dist,-7.0));			
			//printf("i=%d,j=%d, dist ij = %f, force = %f\n",i,j,dist, totalForce);
			force[i*DIM + 0] += totalForce*dx/dist;
			force[i*DIM + 1] += totalForce*dy/dist;
			force[j*DIM + 0] -= totalForce*dx/dist;
			force[j*DIM + 1] -= totalForce*dy/dist;
		}
		if(isNoseHover){				//对Nose-Hover热库
			//eta = 1.0;
			force[i*DIM + 0] += -eta*velocity[i*DIM + 0];
			force[i*DIM + 1] += -eta*velocity[i*DIM + 1];
			//if(eta > 50 ) printf("step=%d,time=%f,fx=%f,fy=%f,x=%6f,y=%6f,vx=%f,vy=%f,eta=%f\n",
			//	step,time,force[i*DIM+0],force[i*DIM+1],pos[i*DIM],pos[i*DIM+1],velocity[i*DIM+0],velocity[i*DIM+1],eta);
		}
		//if(abs(force[i*DIM + 0]>max_force)) max_force = abs(force[i*DIM + 0]);
		//if(abs(force[i*DIM + 1]>max_force)) max_force = abs(force[i*DIM + 1]);
	}
	
	for(int i=0; i<N; i++){
		//printf("i=%d,fx=%f,fy=%f,vx=%f,vy=%f\n",i,force[i*DIM+0],force[i*DIM+1],velocity[i*DIM+0],velocity[i*DIM+1]);
	}
}

// Verlet 算法
void compute_pos(double step = 1e-5){
	for(int i=0;i<N*DIM;i++){
		pre_pos[i] = 2*pos[i]-pre_pos[i]+step*step*force[i];  // Verlet 算法
		velocity[i] = - (pos[i]-pre_pos[i]) / step;
	}
	swap(pre_pos,pos);
}

void compute_energy(){
	double  dist, dx, dy;
	for(int i=0;i<N;i++){
		kinetic_energy[i] = 0.5*velocity[i*DIM+0]*velocity[i*DIM+0]
							+0.5*velocity[i*DIM+1]*velocity[i*DIM+1];
		potential[i] = 0.5*(pos[i*DIM+0]*pos[i*DIM+0]+pos[i*DIM+1]*pos[i*DIM+1]);
		for(int j=0;j<N;j++){
			if (i==j) continue;
			dx = pos[i*DIM + 0] - pos[j*DIM + 0];
			dy = pos[i*DIM + 1] - pos[j*DIM + 1];
			dist = sqrt(dx*dx+dy*dy);
			potential[i] += 4*(pow(dist,-12)-pow(dist,-6));
		}
	}
}

//打印所有粒子的位置信息
void print_pos(int step = 0){
	for(int i=0; i<N; i++){
		printf("step:%d,i=%d,x=%f,y=%f\n",step,i,pos[i*DIM+0],pos[i*DIM+1]);
	}
	cout << endl;
}

//Nose-Hover 热库
//a = f/m -> a=f/m - v*eta
//此部分计算eta
void compute_eta(double timeStep){
	double total_velocity_square = 0.0;
	for(int i=0; i<N; i++){
		total_velocity_square += velocity[i*DIM+0]*velocity[i*DIM+0]
								+velocity[i*DIM+1]*velocity[i*DIM+1];
	}
	eta += timeStep*(total_velocity_square - N*DIM*Temperature)/Ms;	//此处有疑问，公式中的h等于自由度？
}

int main(){
	pos			= new double[N*DIM];
	pre_pos		= new double[N*DIM];
	velocity	= new double[N*DIM];
	force		= new double[N*DIM];
	potential   = new double[N];
	kinetic_energy = new double[N];
	initialize();
	print_pos(0);
	system("pause");
	ofstream out("result.txt");

	for(step=0;step<MAX_STEP;step++){
		time_step = 1e-5;//1.0/max_force;		//计算应该使用的时间步长（取为固定步长）
		if(isNoseHover){
			compute_eta(time_step);				//计算热库摩擦力
		}
		compute_force();					//计算所有粒子受力		
		compute_pos(time_step);				//向前推进一步
		time+=time_step;					//计算此步物理时间
		if(step%10000 == 0)					//每隔10000步存储一次数据
		{			
			compute_energy();
			out << time << " " ;
			for(int in=0; in<N; in++){
				out << pos[in*DIM +0] << " " << pos[in*DIM +1] << " "
					<< velocity[in*DIM +0] << " " << velocity[in*DIM +1] << " "
					<< kinetic_energy[in] << " " << potential[in] << " ";
			}
			out << endl;
			if(step%1000000 == 0)				//每隔1000000显示一次进度
				print_pos(step);
		}
	}
	out.close();
	return 0;
}

