#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <utility>
#include <random>
#include <limits>
#include <iomanip>
#include <chrono>
#include <vector>

const char input_file[] = "Volume.crd";
const char output_file[] = "V_N.txt";

/* a point object has three coordinates in 3D space */
class Point{
protected:
   double _x;
   double _y;
   double _z;
public:
   Point(){}
   Point(const double& x, const double& y, const double&z)
   : _x(x), _y(y), _z(z) {}
   double get_x() const {return _x;}
   double get_y() const {return _y;}
   double get_z() const {return _z;}
};

class Box{
public:
   double x_min = std::numeric_limits<double>::max();
   double x_max = std::numeric_limits<double>::min();
   double y_min = std::numeric_limits<double>::max();
   double y_max = std::numeric_limits<double>::min();
   double z_min = std::numeric_limits<double>::max();
   double z_max = std::numeric_limits<double>::min();
   Box(){};
   double volume() const;
};

/* member function to return the volume of the current Box object */
double Box::volume () const
{
   return (x_max-x_min) * (y_max - y_min) * (z_max - z_min);
}

/* an Box object has a point object plus a radius */
class Atom : public Point
{
private:
   double _r;
public:
   Atom():Point(){};
   Atom(const double& x, const double& y, const double&z, const double& r)
   : Point(x, y, z), _r(r){}
   Atom(const Point& p) : Point(p) {};
   double get_r() const {return _r;}
   
   friend std::ostream& operator<<(std::ostream&, const Atom& atom);
};

/* overloaded insertion opeartion makes the printing of Atom object easier */
std::ostream& operator<<(std::ostream& s,const Atom& atom)
{
   std::cout << atom._x << "\t" << atom._y << "\t" << atom._z << "\t"
   << atom._r << "\n";
   return s;
}

/* calculate the total surface area of the union of all balls*/
double get_total_surface_area(Atom atoms [], const unsigned& N);

/* open the input file so that it is ready to be read */
void open_file(int& size,std::ifstream& fin);

/* open the input file to make reading operation ready */
void process_file(std::ifstream& fin, Box& box, Atom atoms []);

/* get a random number between [min, max] */
double get_random_double(const double& min, const double& max);

/* get a random point within the box */
Point get_random_point(const Box& box);

/* get a random point from the surface of all balls */
Point get_random_point_prime(const Atom& atom);

/* get the distance between two points */
double get_distance(const Point& p, const Atom& a);

/* f(xi) is 1 if the distance from xi to center of any one of these atoms
   is lesser than its radius */
int calculate_fx(const Point& xi, const Atom atoms[], const unsigned& size);

/* f(xi) prime is 1 if the distance from xi to center of ONLY one of these atoms
is same than its radius */
int calculate_fx_prime(const Point& xi, const Atom atoms[], const unsigned& size);


/* calculate the total volume */
void monte_carlo(const double& V, const unsigned int& n,
                                      const Box& box, const Atom atoms [],
                                      const unsigned& size);
/* calculate the toal surface area */
void monte_carlo_like(const double& SA, const unsigned int& N,
                                     const Atom atoms [],
                                     const unsigned& M);

int main(int argc, const char * argv[]) {
   
   srand(static_cast<unsigned int>(time(nullptr)));
   Box box;
   int size; // the number of rows of data in the input file
   std::ifstream fin;
   open_file(size, fin);
   Atom atoms[size];
   unsigned int N = 12900;

   process_file(fin, box, atoms);
   
   fin.close();
   
   double V = box.volume();
      
   monte_carlo(V, N, box, atoms, size);

//   double SA = get_total_surface_area(atoms, size);
//   monte_carlo_like(SA, N, atoms, size);
}

double get_total_surface_area(Atom atoms [], const unsigned& N)
{
   double SA = 0;
   double pi = 3.14159265358979323846;
   for (int i = 0; i < N; i++)
   {
      double curr_sa = 4*pi*atoms[i].get_r();
      SA += curr_sa;
   }
   return SA;
}

void open_file(int& size, std::ifstream& fin)
{
   fin.open(input_file);
   if(!fin){
      std::cerr << "Unable to process the input file: " << input_file << "\n";
      exit(-1);
   }
   fin >> size;
}

void process_file(std::ifstream& fin, Box& box, Atom atoms [])
{
   double x, y, z, r;
   double r_max = std::numeric_limits<double>::min();
   int i = 0;
   
   while (!fin.eof()) {
      
      fin >> x >> y >> z >> r;

      Atom atom(x, y, z, r);
      
      if(r > r_max) r_max = r;
      
      if (x < box.x_min) box.x_min = x;
      if (x > box.x_max) box.x_max = x;
      
      if (y < box.y_min) box.y_min = y;
      if (y > box.y_max) box.y_max = y;
      
      if (z < box.z_min) box.z_min = z;
      if (z > box.z_max) box.z_max = z;
      
      atoms[i++] = atom;
   }
   box.x_min = box.x_min - r_max;
   box.x_max = box.x_max + r_max;
   box.y_min = box.y_min - r_max;
   box.y_max = box.y_max + r_max;
   box.z_min = box.z_min - r_max;
   box.z_max = box.z_max + r_max;
}

double get_random_double(const double& min, const double& max)
{
   std::random_device rd;  //Will be used to obtain a seed for the random number engine
   std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
   std::uniform_real_distribution<> dis(min, max);
   return dis(gen);
}

Point get_random_point(const Box& box)
{
   double random_x = get_random_double(box.x_min, box.x_max);
   double random_y = get_random_double(box.y_min, box.y_max);
   double random_z = get_random_double(box.z_min, box.z_max);
   
   return Point(random_x, random_y, random_z);
}

Point get_random_point_prime(const Atom& atom)
{
   double x_min = atom.get_x() - atom.get_r();
   double x_max = atom.get_x() + atom.get_r();
   double y_min = atom.get_y() - atom.get_r();
   double y_max = atom.get_y() + atom.get_r();
   double random_x = get_random_double(x_min, x_max);
   double random_y = get_random_double(y_min, y_max);
   double x_square = (random_x - atom.get_x()) * (random_x - atom.get_x());
   double y_square = (random_y - atom.get_y()) * (random_y - atom.get_y());
   double r_square = (atom.get_r()) * (atom.get_r());
   double random_z = sqrt(r_square - x_square - y_square);
   return Point(random_x, random_y, random_z);
}

double get_distance(const Point& p, const Atom& a)
{
   double d_x = p.get_x() - a.get_x();
   double d_y = p.get_y() - a.get_y();
   double d_z = p.get_z() - a.get_z();
   
   return sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
}

int calculate_fx(const Point& xi, const Atom atoms[], const unsigned& size)
{
   for(int i = 0; i < size; i++)
   {
      double distance = get_distance(xi, atoms[i]);
      if(distance < atoms[i].get_r())
         return 1;
   }
   return 0;
}

int calculate_fx_prime(const Point& xi, const Atom atoms[], const unsigned& size)
{
   int counter = 0;
   for(int i = 0; i < size; i++)
   {
      double distance = get_distance(xi, atoms[i]);
      if(distance == atoms[i].get_r())
         return 1;
      if(counter == 2)
         return 0;
   }
   return counter;
}

void monte_carlo(const double& V, const unsigned int& n,
                                      const Box& box, const Atom atoms [],
                                      const unsigned& size)
{
   double s = 0, s2 = 0;
   
   // initialize number of random points N
   unsigned int N = n;
   
   double f = 0, f2 = 0, vol = 0, sd = 0;
   
   // create an ouput file
   std::ofstream fout(output_file);
   
   // get the very begining timepoint
   auto start = std::chrono::high_resolution_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::microseconds>(start - start);

   for(int i = 1 ; i <= N; i ++)
   {
      auto begin = std::chrono::high_resolution_clock::now();

      Point xi = get_random_point(box); // position xi at random inside the  box
      int f_of_xi = calculate_fx(xi, atoms, size); // compute f(xi),
      
      // update the sums
      s = s + f_of_xi;
      s2 = s2 + f_of_xi*f_of_xi;
      
      // compute the means
      f = s/i;
      f2 = s2/i;
      
      vol = V * f; // compute the volume
      
      // compute the standard deviation error
      sd = V * sqrt((f2-f*f)/i);
            
      // get stopping timepoint
      auto stop = std::chrono::high_resolution_clock::now();
      duration += std::chrono::duration_cast<std::chrono::microseconds>(stop - begin);
      
      // add results to the outputfile
      fout << std::fixed << std::setprecision(2) << std::setw(10) << std::left
         << vol << "\t"  << std::setw(10) /* estimated V*/
         << sd <<"\t" << std::setw(10) /* estimated V lower bound*/
         << vol - sd << "\t" << std::setw(10)
         << vol + sd << "\t" << std::setw(10)
         << duration.count() << "\t" << std::setw(10)
         << i  << "\n";
   }
}

void monte_carlo_like(const double& SA, const unsigned int& N,
                                        const Atom atoms [],
                                        const unsigned& M)
{
   double s = 0, s2 = 0;
   double f = 0, f2 = 0, sa = 0, sd = 0;
   
   // get the very begining timepoint
   auto start = std::chrono::high_resolution_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::microseconds>(start - start);

   
   // get N randoms points from the surface of all balls
   std::vector<Point> points;
   for(int i = 1 ; i < N; i++)
   {
      Point xi = get_random_point_prime(atoms[i%M]); // position xi at random at the surface
      points.push_back(xi);
   }

   // find probability that a random point falls on the surface
   // of exactly ONE of these balls
   for(int i = 1 ; i < N; i++)
   {
      auto begin = std::chrono::high_resolution_clock::now();
      
      Point xi = points[i];
      
   
      int f_of_xi = calculate_fx_prime(xi, atoms, M); // compute f(xi),
      
      // update the sums
      s = s + f_of_xi;
      s2 = s2 + f_of_xi*f_of_xi;
      
      // compute the means
      f = s/i;
      f2 = s2/i;
      
      sa = SA * f; // compute the volume
      
      // compute the standard deviation error
      sd = sa * sqrt((f2-f*f)/i);
            
      // get stopping timepoint
      auto stop = std::chrono::high_resolution_clock::now();
      duration += std::chrono::duration_cast<std::chrono::microseconds>(stop - begin);
      
      // add results to the outputfile
      std::cout << std::fixed << std::setprecision(2) << std::setw(10) << std::left
         << sa << "\t"  << std::setw(10) /* estimated V*/
         << sd <<"\t" << std::setw(10) /* estimated V lower bound*/
         << sa - sd << "\t" << std::setw(10)
         << sa + sd << "\t" << std::setw(10)
         << duration.count() << "\t" << std::setw(10)
         << i  << "\n";
   }
}

