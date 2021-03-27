#ifndef HEAP_H
#define HEAP_H

//#include "network.h"
#include <vector>
#include <algorithm>

// Generic class used to construct a heap where the ordering is specified by a
// comparator provided to the constructor
template <class T>
class HEAP{
  std::vector<T> elements;
  bool(*comparator) (T,T);
  
public:
  /* Creates a heap with no nodes that uses the provided comparator
   * for ordering. The comparator should return true iff T1 is less
   * than T2. */
  HEAP(bool (*comp) (T,T)){
    comparator = comp;
    elements = std::vector<T> ();
  }
  
  /* Creates a heap using the given nodes and comparator for ordering. 
   * The comparator should return true iff T1 is less than T2. */
  HEAP(std::vector<T> nodes, bool (*comp) (T,T)) {
    comparator = comp;
    elements = nodes;
    std::make_heap(elements.begin(), elements.end(), comparator);
  }
  
  /* Inserts the provided element into the heap */
  void insert(T newNode){
    elements.push_back(newNode);
    std::push_heap(elements.begin(), elements.end(), comparator);
  }

  /* Removes the "largest" value from the heap and returns it */
  T pop(){
    std::pop_heap(elements.begin(), elements.end(), comparator);
    T most = elements.back();
    elements.pop_back();
    return most;
  }

  /* Returns the number of elements contained in the heap */
  int size(){
    return elements.size();
  }
  
  /*Rearrranges the entire heap to ensure that the proper ordering is preserved */
  void reHeapify(){
      std::make_heap(elements.begin(), elements.end(), comparator);
  }
};

#endif
