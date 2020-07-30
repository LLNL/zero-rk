#include <cstdio>
#include <string>
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

} // namespace zerork

