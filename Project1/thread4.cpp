#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <chrono>

using namespace std;

#define NUM_THREADS 4

float calculate(long tid) {
    float val = 0.;
    for (int j = 0; j < 1e9; j++)
        val += 1.0 * (float)tid;

    return val;
}

void* wait(void* t) {
    int i;
    long tid;

    tid = (long)t;
    cout << "Thread with id : " << tid << "  starting " << endl;
    float val = calculate(tid);
    cout << "Result: " << val << endl;

    cout << "Thread with id : " << tid << "  exiting " << endl;
    pthread_exit(NULL);
}

int main() {
    int rc;
    int i;
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    void* status;

    // Initialize and set thread joinable
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    auto start = chrono::high_resolution_clock::now();

    for (i = 0; i < NUM_THREADS; i++) {
        cout << "main() : creating thread, " << i << endl;
        rc = pthread_create(&threads[i], &attr, wait, (void*)i);
        if (rc) {
            cout << "Error:unable to create thread," << rc << endl;
            exit(-1);
        }
    }

    // free attribute and wait for the other threads
    pthread_attr_destroy(&attr);
    for (i = 0; i < NUM_THREADS; i++) {
        rc = pthread_join(threads[i], &status);
        if (rc) {
            cout << "Error:unable to join," << rc << endl;
            exit(-1);
        }
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Duration (threaded): " << duration.count() << " ms" << endl;

    start = chrono::high_resolution_clock::now();
    cout << calculate(1) << endl;
    cout << calculate(2) << endl;
    cout << calculate(3) << endl;
    cout << calculate(4) << endl;
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Duration (linear): " << duration.count() << " ms" << endl;

    pthread_exit(NULL);
}

