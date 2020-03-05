/* Code to compute alpha value for single connected component
   applied to run on each frame of an XTC file.
   Code by Peter Kasson, based on CGAL documentation.
*/
#include <gflags/gflags.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "Alpha_shape_3mod.h"
#include <xtcio.h>
#include <smalloc.h>
#include <fstream>
#include <list>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;
typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;
typedef Gt::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;

class trajObj {
public:
  std::list<Point> get_frame();
  trajObj(const char *fname);
  bool next();
  ~trajObj();
private:
  std::list<Point> pointlist;
  matrix box;
  rvec *x;
  real time;
  int step;
  int natoms;
  t_fileio *fio;
public:
  real get_time() {return time;}
};

trajObj::trajObj(const char *fname) {
  // Open file and read first frame
  real prec;
  gmx_bool bOK;
  x = NULL;
  natoms = 0;
  fio = open_xtc(fname, "r");
  read_first_xtc(fio, &natoms, &step, &time, box, &x, &prec, &bOK);
  if (!bOK) {
    std::cerr << "Error reading first frame of xtc!";
  }
}

trajObj::~trajObj() {
  // Destructor
  if (x) {
    sfree(x);
  }
  close_xtc(fio);
}

std::list<Point> trajObj::get_frame() {
  // Get the frame as a list of Points
  pointlist.clear();
  for (int i = 0; i < natoms; i ++) {
    pointlist.push_back(Point(x[i][0], x[i][1], x[i][2]));
  }
  return pointlist;
}

bool trajObj::next() {
  // Read next frame
  gmx_bool bOK;
  real prec;
  if (read_next_xtc(fio, natoms, &step, &time, box, x, &prec, &bOK)) {
    return true;
  } else {
    return false;
  }
}

DEFINE_bool(enumerate, false, "Enumerate component sizes");

float alpha_threshold(std::list<Point> lp, float alpha_val,
                      int thresh=50, bool debug=false)
{
  //DEBUG
  if (debug) {
    std::cerr << "Number of points: " << lp.size() << std::endl;
  }
  // compute alpha shape
  Alpha_shape_3 as(lp.begin(),lp.end(), alpha_val);
  std::list<int> components = as.count_solid_components(alpha_val);
  //threshold components
  int component_ctr=0;
  for (std::list<int>::iterator it = components.begin(); it!=components.end(); it++)
  {
    if (*it >= thresh) {
      component_ctr++;
    }
    if (FLAGS_enumerate) {
     std::cout << *it << '\t';
    }
  }

  // find optimal alpha value
  // Alpha_iterator opt = as.find_optimal_alpha(1);
  // as.set_alpha(*opt);
  // assert(as.number_of_solid_components() == 1);
  return component_ctr;
}

DEFINE_string(xtcfile, "", "Input file.");
DEFINE_double(alpha, 0.12, "Alpha cutoff");
DEFINE_int32(thresh, 50, "Minimum simplices in a component");

int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  trajObj traj(FLAGS_xtcfile.c_str());
  // Ultimately better output here
  std::cout << "Time\tComponents" << std::endl;
  do {
    std::cout << traj.get_time() << '\t'
              << alpha_threshold(traj.get_frame(),
                                 static_cast<float>(FLAGS_alpha), FLAGS_thresh)
              << std::endl;
  } while (traj.next());
  return 0;
}

