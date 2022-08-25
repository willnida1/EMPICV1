/*
  Three-vector class including mathematical operations and IO
 */

//#include <cmath>
#include <iostream>
#include <fstream>

class ThreeVec
{
private:
  double coord_[3]; // Private data members e.g. x,y,z

public:
  // Default constructor
  ThreeVec()
  {
    for (int i = 0; i < 3; ++i)
      coord_[i] = 0.0;
  }

  // Cartesian constructor
  ThreeVec(double x, double y, double z)
  {
    coord_[0] = x;
    coord_[1] = y;
    coord_[2] = z;
  }
  // Access function for x coordinate
  double getX()
  {
    return coord_[0];
  }

  // Access function for y coordinate
  double getY()
  {
    return coord_[1];
  }

  // Access function for z coordinate
  double getZ()
  {
    return coord_[2];
  }
  // Access function for ith coordinate
  double get(int i)
  {
    return coord_[i];
  }

  // Modifier method for x coordinate
  void setX(double value)
  {
    coord_[0] = value;
  }

  // Modifier method for y coordinate
  void setY(double value)
  {
    coord_[1] = value;
  }

  // Modifier method for z coordinate
  void setZ(double value)
  {
    coord_[2] = value;
  }
  // Modifier method for ith coordinate -> INSERT
  void set(int i, double value)
  {
    coord_[i] = value;
  }
  // Alternative modifier method for ith coordinate -> ADD
  void inc(int i, double value)
  {
    coord_[i] += value;
  }

  // Square the threevector
  double square()
  {
    double answer = 0.0;
    for (int i = 0; i < 3; ++i)
      answer += coord_[i] * coord_[i];
    return answer;
  }

  // Magnitude of the threevector
  double mag()
  {
    return sqrt(square());
  }

  /*
    Overload the operators +,-,* and ^ to represent vector operations
  */
  // Addition
  ThreeVec operator+(ThreeVec vec)
  {
    ThreeVec ans(vec.getX() + coord_[0],
                 vec.getY() + coord_[1],
                 vec.getZ() + coord_[2]);
    return ans;
  }

  // Subtraction
  ThreeVec operator-(ThreeVec vec)
  {
    ThreeVec ans(coord_[0] - vec.getX(),
                 coord_[1] - vec.getY(),
                 coord_[2] - vec.getZ());
    return ans;
  }

  // Increment
  ThreeVec operator+=(ThreeVec vec)
  {
    coord_[0] += vec.getX();
    coord_[1] += vec.getY();
    coord_[2] += vec.getZ();
    return *(this); // not necessary but returns modified ThreeVec
  }

  // Scalar multiplication
  ThreeVec operator*(double value)
  {
    ThreeVec ans(coord_[0] * value,
                 coord_[1] * value,
                 coord_[2] * value);
    return ans;
  }
  // Scalar division
  ThreeVec operator/(double value)
  {
    ThreeVec ans(coord_[0] / value,
                 coord_[1] / value,
                 coord_[2] / value);
    return ans;
  }
  // Vector multiplication - dot product
  double operator*(ThreeVec vec)
  {
    double ans = 0.0;
    for (int i = 0; i < 3; ++i)
      ans += coord_[i] * vec.get(i);
    return ans;
  }

  // Vector multiplication - cross product
  ThreeVec operator^(ThreeVec vec)
  {
    ThreeVec ans(coord_[1] * vec.getZ() - coord_[2] * vec.getY(),
                 coord_[2] * vec.getX() - coord_[0] * vec.getZ(),
                 coord_[0] * vec.getY() - coord_[1] * vec.getX());
    return ans;
  }

  // Print to screen
  void print()
  {
    for (int i = 0; i < 3; ++i)
      std::cout << coord_[i] << '\t';
    std::cout << '\n';
  }

  // Print to file
  void print(std::ofstream &fout)
  {
    for (int i = 0; i < 3; ++i)
      fout << coord_[i] << '\t';
    fout << '\n';
  }
};
