#include<ctime>
#include<iostream>
using std::cout;
using std::endl;
/*测试排序算法运行时间*/
void testSort(void(*sort)(int[], int), int arr[], int length)
{
    clock_t startTime = clock();
    sort(arr, length);
    clock_t endTime = clock();
    cout << double(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}
/*插入排序*/
void insertionSort(int arr[], int n)
{
    for (int i = 1; i < n; i++)
    {
        int e = arr[i];
        int j;
        for (j = i; j >0 && e < arr[j - 1]; j--)
            arr[j] = arr[j - 1];
        arr[j] = e;
    }

}

int main()
{
    int arrSize = 10000;
    int *arr = new int[arrSize];
    srand(time(NULL));
    for (int i = 0; i < arrSize; i++)
    {
        arr[i] = rand() % arrSize;
    }
    testSort(insertionSort, arr, arrSize);

    system("pause");
    return 0;
} 