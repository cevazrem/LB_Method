#include <iostream>
#include <vector>
#include <fstream>

//const int m = 100; // slab length&
//const int time_step = 200;

struct all_data{
    double f0, f1, f2;
    double rho, feq, x;
};

class lb_method {
    private:
        std::vector <all_data> data;
        double dt, dx, csq, alpha, omega, left_temp, start_temp; 
        int m, time_step; 
		std::ofstream result;
    public:
    
    lb_method() {
		m = 100;
		time_step = 200;
		data.resize(m);
		result.open("Result.txt");
        dt = 1;
        dx = 1;
        data[0].x = start_temp;
        for (int i = 1; i < m; ++i) {
            data[i].x = data[i-1].x + dx;
        }
        csq = (dx * dx)/(dt * dt);
        alpha = 0.25;
        omega = 1 / (alpha / (dt * csq) + 0.5);
        left_temp = 1;
        start_temp = 0;
        for (int i = 0; i < m ; ++i) {
			data[i].rho = 0;
			data[i].f1 = 0.5 * data[i].rho;
			data[i].f2 = 0.5 * data[i].rho;
		}
    }
    
    lb_method(double length, double cur_time, double _alpha, double _left_temp, double _start_temp, int _m, int _time_step) {
		m = _m;
		time_step = _time_step;
		data.resize(m);
		result.open("Result.txt");
        dt = cur_time/time_step;
        dx = length/m;
        data[0].x = start_temp;
        for (int i = 1; i < m; ++i) {
            data[i].x = data[i-1].x + dx;
        }
        csq = (dx * dx)/(dt * dt);
        alpha = _alpha;
        omega = 1 / (alpha / (dt * csq) + 0.5);
        left_temp = _left_temp;
        start_temp = _start_temp;
        for (int i = 0; i < m ; ++i) {
			data[i].rho = 0;
			data[i].f1 = 0.5 * data[i].rho;
			data[i].f2 = 0.5 * data[i].rho;
		}
    }
    
    void collision() {
        for (int i = 0; i < m; ++i) {
            data[i].rho = data[i].f1 + data[i].f2;
            data[i].feq = 0.5 * data[i].rho;
            data[i].f1 = (1 - omega) * data[i].f1 + omega * data[i].feq;
            data[i].f2 = (1 - omega) * data[i].f2 + omega * data[i].feq;
        }
    }
    
    void streaming() {
		for (int i = 1; i < m - 1; ++i) {
            data[m - i - 1].f1 = data[m - i - 2].f1;
            data[i - 1].f2 = data[i].f2;
        }
	}
	
	void bound_cond() {
		data[0].f1 = left_temp - data[0].f2;
        data[m - 1].f1 = data[m - 2].f1;
        data[m - 1].f2 = data[m - 2].f2;
	}
	
	void write_file() {
		for (int i = 0 ; i < m; ++i) {
			result << data[i].x << ' ' << data[i].rho << std::endl;
		}
	}
	
	void calculate() {
		for (int i = 0; i < time_step; ++i) {
			collision();
			streaming();
			bound_cond();
		}
	}
};


/*
int main() {
    int i;
    int m = 100;
    std::ofstream result("Result.txt");
    if (!result.is_open()) {
        std::cout << "Cant open the result file";
        return 1;
    }
    std::vector <all_data> data(m); 
    double dt = 1;
    double dx = 1;
    data[0].x = 0;
    for (int i = 1; i < m; ++i) {
        data[i].x = data[i-1].x + dx;
    }
    double csq = (dx * dx)/(dt * dt);
    double alpha = 0.25;
    double omega = 1 / (alpha / (dt * csq) + 0.5);
    int time_stap = 200;
    double left_temp = 1;
    double start_temp = 0;
    
    for (int i = 0; i < m ; ++i) {
        data[i].rho = 0;
        data[i].f1 = 0.5 * data[i].rho;
        data[i].f2 = 0.5 * data[i].rho;
    }    

    for (int tt = 1; tt < time_stap; ++tt) {
        for (int i = 0; i < m; ++i) {
            data[i].rho = data[i].f1 + data[i].f2;
            data[i].feq = 0.5 * data[i].rho;
            data[i].f1 = (1 - omega) * data[i].f1 + omega * data[i].feq;
            data[i].f2 = (1 - omega) * data[i].f2 + omega * data[i].feq;
        }
        
        for (int i = 1; i < m - 1; ++i) {
            data[m - i - 1].f1 = data[m - i - 2].f1;
            data[i - 1].f2 = data[i].f2;
        }

        data[0].f1 = left_temp - data[0].f2;
        data[m - 1].f1 = data[m - 2].f1;
        data[m - 1].f2 = data[m - 2].f2;
    }
    
    for (int i = 0 ; i < m; ++i) {
        result << data[i].x << ' ' << data[i].rho << std::endl;
    }
    return 0;
}
*/

int main(int argc, char **argv) {
	int mode;
	std::cout << "Enter mode. 1-model, 2-config" << std::endl;
	std::cin >> mode;
	lb_method solution;
	if (mode == 1)
		solution = lb_method();
	else {
		double length, cur_time, alpha, left_temp, start_temp;
		int m, time_step;
		std::cout << "Enter lenght of layer:";
		std::cin >> length;
		std::cout << std::endl << "Enter time point:";
		std::cin >> cur_time;
		std::cout << std::endl << "Enter alpha:";
		std::cin >> alpha;
		std::cout << std::endl << "Enter left side temperature:";
		std::cin >> left_temp;
		std::cout << std::endl << "Enter system temperature t=0:";
		std::cin >> start_temp;
		std::cout << std::endl << "Enter ammount of x lattice:";
		std::cin >> m;
		std::cout << std::endl << "Enter ammount of time steps:";
		std::cin >> time_step;
		std::cout << std::endl;
		solution = lb_method(length, cur_time, alpha, left_temp, start_temp, m, time_step);
		
	}
	solution.calculate();
	solution.write_file();
	return 0;
}
