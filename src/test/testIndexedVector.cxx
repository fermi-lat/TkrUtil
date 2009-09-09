// IndexedVector.cxx : Defines the entry point for the console application.
//

#include <iostream>
#include <stdexcept>
#include "TkrUtil/IndexedVector.h"
using namespace std;

int main(int /* argc*/, char* /*argv[]*/)
{
    typedef IndexedVector<int> intVec;
    intVec myVec;

    myVec.setDims(10,10,10,10,10);

    cout << "Test that correct values are filled and retrieved" << endl;
    cout << "We expect 100,000 correct trials" << endl;
    int i, j, k, m, n;
    for(i=0;i<10;++i) {
        for(j=0;j<10;++j) {
            for (k=0;k<10;++k) {
                for (m=0;m<10;++m) {
                    for (n=0;n<10;++n) {
                        myVec(i,j,k,m,n) = (i+1)*(j+1)*(k+1)*(m+1)*(n+1);
                    }
                }
            }
        }
    }
    
    int okCount = 0;
    for(i=0;i<10;++i) {
        for(j=0;j<10;++j) {
            for (k=0;k<10;++k) {
                for (m=0;m<10;++m) {
                    for (n=0;n<10;++n) {
                        int ind = myVec.getIndex(i, j, k, m, n);
                        int test = (i+1)*(j+1)*(k+1)*(m+1)*(n+1);
                        if (myVec[ind]!=test || myVec[ind]!=myVec(i,j,k,m,n)) {
                            cout << "oopsla!" << endl;
                        } else {
                            okCount++;
                        }
                    }
                }
            }
        }
    }
   cout << "Done, " << okCount << " elements retrieved correctly" << endl;

    cout << "Now two calls with invalid arguments, exceptions should be thrown" << endl;
    try {
        cout << myVec(-1,0,0,0,0) << endl;
    } catch (invalid_argument e) {
        cout << "exception caught: " << e.what() << endl;
    }

     try {
        cout << myVec(0,0,0,10,0) << endl;
    } catch (invalid_argument e) {
        cout << "exception caught: " << e.what() <<  endl;
    }
    return true;
}

