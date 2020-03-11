#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <utility>
#include <random>
#include <limits>

const char input_file[] = "Volume.crd";
const char output_file[] = "VandN.txt";

using namespace std;

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

double Box::volume () const
{
   return (x_max-x_min) * (y_max - y_min) * (z_max - z_min);
}

class Atom : public Point
{
private:
   double _r;
public:
   Atom(){};
   Atom(const double& x, const double& y, const double&z, const double& r)
   : Point(x, y, z), _r(r) {}
   bool operator>(const Atom&) const;
   bool operator<(const Atom&) const;
   bool operator>=(const Atom&) const;
   bool operator<=(const Atom&) const;
   double get_r() const {return _r;}
   
   friend ostream& operator<<(ostream&, const Atom& atom);
};

/* enable > comparasion between two atoms */
bool Atom::operator>(const Atom& atom) const {
   if(_x > atom._x && _y > atom._y && _z > atom._z)
      return true;
   else
      return false;
}

/* enable < comparasion between two atoms */
bool Atom::operator<(const Atom& atom) const {
   return atom > *this;
}

/* enable >= comparasion between two atoms */
bool Atom::operator>=(const Atom& atom) const {
   if(_x >= atom._x && _y >= atom._y && _z >= atom._z)
      return true;
   else
      return false;
}

/* enable <= comparasion between two atoms */
bool Atom::operator<=(const Atom& atom) const {
   return atom >= *this;
}

ostream& operator<<(ostream& s,const Atom& atom)
{
   cout << atom._x << "\t" << atom._y << "\t" << atom._z << "\t"
   << atom._r << "\n";
   return s;
}

ifstream open_file(int& size);

/* open the input file to make reading operation ready */
void process_file(ifstream& fin, Box& box, Atom atoms []);

/* get a random number between [min, max] */
double get_random_double(const double& min, const double& max);

/* get a random point within the box */
Point get_random_point(const Box& box);

/* get the distance between two points */
double get_distance(const Point& p, const Atom& a);

/* f(xi) is 1 if the distance from xi to center of any atom inside the box is
 lesser than its radius */
int calculate_fx(const Point& xi, const Atom atoms[], const unsigned& size);

/* calculate the total volume */
std::pair <double,double> monte_carlo(const double& V, const unsigned int& n,
                                      const Box& box, const Atom atoms [],
                                      const unsigned& size);

int main(int argc, const char * argv[]) {
   
   srand(time(NULL));
   Box box;
   int size; // the number of rows of data in the input file
   ifstream fin = open_file(size);
   Atom atoms[size];
   
   process_file(fin, box, atoms);
   
   // create an ouput file
   ofstream fout(output_file);
   
   double V = box.volume();
   unsigned int N = 100;
   
   //   for (int i = 0 ; i < size; i++) {
   //      cout << atoms[i].get_x() << "\t" << atoms[i].get_y() << "\t" << atoms[i].get_z() << endl;
   //   }
   
   for (int i = 1; i <= N; i++) {
      std::pair <double,double> results = monte_carlo(V,N, box, atoms, size);
      fout << fixed << setprecision(2) << setw(10);
      fout << results.first << " += " << results.second << "\t" << i << "\n";
   }
}

ifstream open_file(int& size)
{
   ifstream fin(input_file);
   if(!fin){
      cerr << "Unable to process the input file: " << input_file << "\n";
      exit(-1);
   }
   fin >> size;
   return fin;
}

void process_file(ifstream& fin, Box& box, Atom atoms [])
{
   double x, y, z, r;
   double r_max = std::numeric_limits<double>::min();
   int i = 0;
   
   while (!fin.eof()) {
      
      fin >> x >> y >> z >> r;
      //      cout << x << "\t"  << y << "\t"  << z << "\t"  << r << "\n";
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

std::pair <double,double> monte_carlo(const double& V, const unsigned int& n,
                                      const Box& box, const Atom atoms [],
                                      const unsigned& size)
{
   double s = 0, s2 = 0;
   
   // initialize number of random points N
   unsigned int N = n;
   
   for(int i = 0 ; i <= N; i ++)
   {
      // position xi at random inside the  box
      Point xi = get_random_point(box);
      
      // compute f(xi),
      int f_of_xi = calculate_fx(xi, atoms, size);
      
      // update the sums
      s = s + f_of_xi;
      s2 = s2 + f_of_xi*f_of_xi;
   }
   
   // comput the means
   double f = s/N;
   double f2 = s2/N;
   
   // compute the volume
   double vol = V * f;
   
   // compute the standard deviation error estimate
   double sd = V * sqrt((f2-f*f)/N);
   
   return make_pair(vol, sd);
}
