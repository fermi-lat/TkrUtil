/// @file IndexedVector.h

/** @class IndexedVector
@brief "Multi-dimensional" vector based on std::vector (should be in some more general package)

@author Leon Rochester
$Header$
*/

#include <vector>

template <class T> class IndexedVector : public std::vector<T>
{
    /* 
    class to allow us to generate effective multi-dimensional vectors
    Data is stored in a simple vector, but instantiated as, for example:
    IndexedVector<float> myArray(4,3,2);
    and accessed as:
    myArray(i, j, k) = 100.f;
    or
    float x = myArray(2,0,1);

    Hand-coded for NUMDIMS possible dimensions (currently 5)

    throws an invalid_argument exception at instantiation or access

    to do:
    * specialize for vectors of pointers, to handle possible delete objects pointed to??
    or do we let the user do this in the destructor of her class?
    */

public:

    enum { NUMDIMS = 5 };

    /// instantiates, but doesn't construct the array
    IndexedVector(): m_nDim(0), m_nElements(0), m_valid(false), m_checkRange(true) {}

    /// instantiates and constructs the array
    IndexedVector(int dim0, int dim1=0, int dim2=0, int dim3=0, int dim4=0)
    {
        if (!setDims(dim0, dim1, dim2, dim3, dim4)) {
            throw std::invalid_argument("IndexedVector: instantiation");
        }
    }

    /// constructs the array after default instantiation
    // returns bool status: true = valid args, false = invalid args
    bool setDims(int dim0, int dim1=0, int dim2=0, int dim3=0, int dim4=0) {

        m_dim[0] = dim0;
        m_dim[1] = dim1;
        m_dim[2] = dim2;
        m_dim[3] = dim3;
        m_dim[4] = dim4;

        m_valid = true;

        // vector is invalid if first dimension is zero
        //   or if there is a negative dimension before the first zero 
        if (m_dim[0]<=0) {
            m_valid = false;
        } else {
            m_nElements = m_dim[0];
            m_nDim = 1;
            for (int i=1;i<NUMDIMS;++i) {
                if (m_dim[i]<0) {
                    m_valid = false;
                    break;
                }
                if (m_dim[i]==0) {
                    break;
                }
                m_nElements *= m_dim[i];
                m_nDim++;
            }
        }
        if (!m_valid) {
            m_nDim = 0;
            m_nElements = 0;
        }
        (*this).resize(m_nElements);
        return m_valid;
    }

    void setValue(T val) {this->assign(m_nElements, val); }

    /// returns a ref to an element by index into each dimension
    T& operator ()(int i, int j=0, int k=0, int m=0, int n=0) {
        return (*this)[getIndex(i, j, k, m, n)];
    }

    bool checkIndex(int i, int j=0, int k=0, int m=0, int n=0) const {
        /// checks the range of the arguments
        // note: "break"s are conditinal on purpose; code sweeps through
        //     currently defined dimensions
        bool valid = false;
        switch (m_nDim) {
            case 5:
                if (n<0 || n>=m_dim[4]) break;
            case 4:
                if (m<0 || m>=m_dim[3]) break;
            case 3:
                if (k<0 || k>=m_dim[2]) break;
            case 2:
                if (j<0 || j>=m_dim[1]) break;
            case 1:
                valid = (i>=0 && i<m_dim[0]);
                break;
            default:
                valid = false;
        }
        return valid;
    }

    /// returns the overall index given the individual indices
    //  if range checking is on, throws an exception for invalid args
    //  since IndexedVector is an std::vector, elements can be accessed like:
    //        int index = myArray.getIndex(i, j, k);
    //        myArray[index] = 0;

    int getIndex(int i, int j=0, int k=0, int m=0, int n=0) const {
        m_validIndex = false;

        // range checking...
        if (m_checkRange) {
           if(!checkIndex(i, j, k, m, n)) throw std::invalid_argument("IndexedVector: access");
        }

        // this could be a loop, but it would probably be slower
        switch(m_nDim){
            case 1:  
                return i;
            case 2: 
                return i*m_dim[1]+j;
            case 3:  
                return (i*m_dim[1]+j)*m_dim[2]+k;
            case 4:  
                return ((i*m_dim[1]+j)*m_dim[2]+k)*m_dim[3]+m ;
            case 5:  
                return (((i*m_dim[1]+j)*m_dim[2]+k)*m_dim[3]+m)*m_dim[4]+n;
            default: // needn't happen if the status is checked when dims are set
                throw std::invalid_argument("IndexedVector: access");
        }
    }

    /// sets range-checking off or on; default is on
    void setRangeCheck(bool check) { m_checkRange = check;}

    /// check for valid vector; false means illegal dimension
    bool isValid() const           { return m_valid; }
    // check for valid index; false means no checking done or index out of range
    bool isValidIndex() const      { return m_validIndex; }
    /// quick way to return the number of elements
    int getNumElements() const     { return m_nElements; }
    /// get the dimensions
    int getDim(int i) const        { return m_dim[i]; }

private:

    /// array to store dimensions
    int          m_dim[NUMDIMS];
    /// current number of dimensions set
    int          m_nDim;
    /// total number of elements
    int          m_nElements;
    /// vector valid flag
    bool         m_valid;
    /// index valid flag
    mutable bool m_validIndex;
    /// check range flag
    bool         m_checkRange;
};

