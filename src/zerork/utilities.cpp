#include <cstdio>
#include <string>

#ifdef _WIN32
#include <windows.h>
#include <profileapi.h>
#else
#include <sys/time.h>
#endif

#include "utilities.h"

namespace zerork {

void upperToLower(std::string s)
{
  int len=s.length();
  int j;
  char c[]=" ";
  char lowerCaseDiff='a'-'A';

  for(j=0; j<len; j++)
 {
      c[0]=s[j];
      if(isupper(c[0])!=0)
      {
          c[0]+=lowerCaseDiff;
          s.replace(j,j+1,c);
        }
  }
  std::printf("upperToLower() = %s|\n",s.c_str());
}

void upperToLower(char *s)
{
  int j=0;
  const char lowerCaseDiff='a'-'A';

  if(s!=NULL)
  {
      while(s[j]!='\0')
      {
          if(('A' <= s[j]) && (s[j] <= 'Z'))
            {s[j]+=lowerCaseDiff;}
          j++;
      }
  }
}

double getHighResolutionTime(void)
{
#ifndef _WIN32
    struct timeval tod;

    gettimeofday(&tod, NULL);
    double time_seconds = (double) tod.tv_sec + ((double) tod.tv_usec / 1000000.0);
#else
    static LARGE_INTEGER Frequency;
    Frequency.QuadPart = 0;
    if(Frequency.QuadPart == 0) {
        QueryPerformanceFrequency(&Frequency);
    }

    LARGE_INTEGER Counts;
    QueryPerformanceCounter(&Counts);

    double time_seconds = ((double)Counts.QuadPart) / Frequency.QuadPart;
#endif
    return time_seconds;
}

} // namespace zerork

