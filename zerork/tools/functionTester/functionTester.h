#ifndef ZERORK_TOOLS_FUNCTIONTESTER_H
#define ZERORK_TOOLS_FUNCTIONTESTER_H

// TODO replace GetLine with C (stdio.h) getline() function
#define MAX_LINE_LEN 1024
int GetLine(FILE *InFile, char *ReadLine, char UntilChar, int MaxChar);

#endif
