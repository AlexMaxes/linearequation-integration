// Timer.h
//
////////////////////////////////////////////////////////////

#ifndef _QT_TIMER_H_INCLUDED_
#define _QT_TIMER_H_INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
///////////////////////////////////////////////////////////
#include <windows.h>
//for LARGE_INTEGER, QueryPerformanceCounter()
//and QueryPerformanceFrequency()
///////////////////////////////////////////////////////////

class Timer  
{
private:
  LARGE_INTEGER begin_;
  LARGE_INTEGER counter_;
  LARGE_INTEGER freq_;  //为准确考虑还是保存频率. 调用seconds()省时间
  bool running_;

public:
  Timer(void): running_(false)
  {
    QueryPerformanceFrequency(&freq_);
    reset();
  }
  void reset(void)
  {
    begin_.QuadPart =  counter_.QuadPart = 0;
	running_ = false;
  }
  void start(void)
  {
    if(running_ == false)
    { 
      running_ = true;
      ::QueryPerformanceCounter(&begin_);
    }
  }
  void stop(void)
  {
    LARGE_INTEGER end;
    if(running_ == true)
    {
      ::QueryPerformanceCounter(&end);
      counter_.QuadPart += end.QuadPart - begin_.QuadPart;
      //begin_.QuadPart = 0;
      running_ = false;
    }
  }
  double seconds(void) const
  {
    if(running_ == false)
      return counter_.QuadPart / static_cast<double>(freq_.QuadPart);
      
    LARGE_INTEGER end;
    ::QueryPerformanceCounter(&end);
    return (counter_.QuadPart + (end.QuadPart - begin_.QuadPart))
         / static_cast<double>(freq_.QuadPart);
  }
  bool running(void) const
  {
    return running_; 
  }
  unsigned long frequency(void) const
  {
    return freq_.LowPart;
  }
  void clear(void)
  {
    reset(); 
  } 
  void go(void)
  {
    start();
  }
  void pause(void) 
  {
    stop(); 
  } 
};//~class Timer

#endif //~ #infdef _QT_TIMER_H_INCLUDED_
//~Timer.h
