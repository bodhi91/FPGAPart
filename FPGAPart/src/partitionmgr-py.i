//
%{
#include "par/PartitionMgr.h"
%}

%include "../../Exception-py.i"

%include <std_vector.i>
%include <std_string.i>

namespace std
{
  %template(vector_int) std::vector<int>;
  %template(vector_float) std::vector<float>;
}

%include "par/PartitionMgr.h"
