#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	//printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	//clprintf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}



// implement heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapify(int arr[], int n, int i) {
    int largest = i; // Initialize largest element as the root 
    int left = 2 * i + 1; // each parent at the left has max of 2 children
    int right = 2 * i + 2; 

    //checks if the left child of root is larger than the current node
    if (left < n && arr[left] > arr[largest])
        largest = left;//left child is now the biggest

    //checks if the right child of root is larger than the current node
    if (right < n && arr[right] > arr[largest])
        largest = right;//right child is now the biggest

    // If largest is not root
    if (largest != i) {
        // Swap to get root to hold the largest element again
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;

      
        heapify(arr, n, largest);
    }
}

void heapSort(int arr[], int n) {
     // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to end
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    } 
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated

void merge(int pData[], int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    // Alloc memory for temporary arrays
    int *L = (int *)Alloc(n1 * sizeof(int));
    int *R = (int *)Alloc(n2 * sizeof(int));

    // Copy data to temporary arrays L[] and R[]
    for (i = 0; i < n1; i++)
        L[i] = pData[l + i];
    for (j = 0; j < n2; j++)
        R[j] = pData[m + 1 + j];

    // Merge the temp arrays back into pData[l..r]
    i = 0; j = 0; k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            pData[k] = L[i];
            i++;
        } else {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of L[], if any
    while (i < n1) {
        pData[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of R[], if any
    while (j < n2) {
        pData[k] = R[j];
        j++;
        k++;
    }

    // Free the temporary arrays
    DeAlloc(L);
    DeAlloc(R);
}

// Main function to implement mergeSort
void mergeSort(int pData[], int l, int r) {
    if (l < r) {
        // Find the middle point to divide the array into two halves
        int m = l + (r - l) / 2;

        // Call mergeSort on the first half
        mergeSort(pData, l, m);

        // Call mergeSort on the second half
        mergeSort(pData, m + 1, r);

        // Merge the sorted halves
        merge(pData, l, m, r);
    }
}


// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{
	int i, item, j;
    for (i = 1; i < n; i++)
    {
         item = pData[i];

        /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
          for(j=i-1; j>=0; j--)
          {
              if(pData[j]>item)
                pData[j+1] = pData[j];
              else
                break;

          }
          pData[j+1] = item;
    }

}

//implement bubble sort
//extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
	int temp;
    for (int i = 0; i < n - 1; i++) {
        // Last i elements are already in place
        for (int j = 0; j < n - i - 1; j++) {
            if (pData[j] > pData[j + 1]) {
                // swap pData[j] and pData[j+1]
                temp = pData[j];
                pData[j] = pData[j + 1];
                pData[j + 1] = temp;
            }
        }
    }
}

//implement selection sort
//extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
	int i, j, min;

    for (i = 0; i < n - 1; i++) {
        min = i;
        for (j = i + 1; j < n; j++) {
            if (pData[j] < pData[min]) {
                min = j;
            }
        }

        if (min != i) {
            int temp = pData[min];
            pData[min] = pData[i];
            pData[i] = temp;

            
        }
		
    }
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData) {
    FILE* inFile = fopen(inputFileName, "r");
    int dataSz = 0;
    

    if (inFile) {
        fscanf(inFile, "%d\n", &dataSz);
        *ppData = (int *)Alloc(sizeof(int) * dataSz);
        if (*ppData == NULL) {
            printf("Cannot allocate memory\n");
            fclose(inFile);
            return -1;
        }

        for (int i = 0; i < dataSz; ++i) {
            fscanf(inFile, "%d", &(*ppData)[i]);
        }
          
        fclose(inFile);
    } else {
        printf("Could not open file\n");
        return -1;
    }

    return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("Merge Sort:\n");
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

        printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}