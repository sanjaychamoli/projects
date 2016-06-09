#include "Solver.h"


int main ( int argc, char **argv ) {

    size_t nx   = 50;    size_t ny   = 50;    size_t c    = 10000; size_t threads = 32;

    if ( argc != 4 ) {  std::cout << " Usage: ./rbgs nx ny c " << std::endl; abort();}
    else {  nx = std::atoi(argv[1]);    ny = std::atoi(argv[2]);    c = std::atoi(argv[3]); }

    std::ofstream file;
    file.open(std::to_string(nx) + "_timing.dat");
    file << "Size: ," << nx << "," << ny;
    file << "\nThreads,Time\n";

    std::cout << "========================\n";
    std::cout << "Size: " << nx << " " << ny << std::endl;
    std::cout << "------------------------\n";
    for(size_t t = 1, i = 0; t <= threads; t*= 2, ++i)
    {
    vec_t u( (nx+1)*(ny+1), 0.0 ), f( (nx+1)*(ny+1), 0.0 ), r( (nx+1)*(ny+1), 0.0 );

    db(u, nx, ny);
	rhs(f, nx, ny);
    vec_t time_data;
    siwir::Timer time;
        time.reset();
            std::cout << std::endl <<"Residual: " << rbgs(u, f, r, nx, ny, c, t);
            time_data.push_back(time.elapsed());
            std::cout << std::endl << "Time : " << time_data.back() << std::endl;
            std::cout << "------------------------\n";
            file << t << "," << time_data.back() << "\n";

            if (t == threads)
            {
                writetofile(u, nx, ny, "solution.txt");
            }
    }
    file.close();
    return 0;
}
