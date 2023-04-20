#include <iostream>
#include <vector>
#include <string>

using namespace std;

float* euler(int N, float y_0, float t_end)
{
    float del_t = t_end/N;
    float array[N];
    for(int n = 0; n<N; n++)
    {
        y_0 = y_0*(1-del_t);
        array[n] = y_0;
    }
    return array;
}

int main()
{
    int N = 100;
    float *test;
    test = euler(N,1,10);
    cout << test;





return 0;
}