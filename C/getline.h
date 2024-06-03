#ifndef H_GETLINE
#define H_GETLINE

extern int getline_reserve (char **lineptr, size_t *n, FILE *stream);
extern int getdelim_reserve (char **lineptr, size_t *n, int delimiter, FILE *fp);

#endif
