#include <stdio.h>
#include <stdlib.h>

int main ()
{
    int i = 0;
    int k = 0;
    int *j;

    j = &k;

    if(j == &i){
        printf("%ls", j);
    }

    printf("%d", &j);
}