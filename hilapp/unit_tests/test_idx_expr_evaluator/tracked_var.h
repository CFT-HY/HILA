
#include <set>

template <typename T>
class TrackedVar {
  public:
    T val;
    std::set<T> values;

    operator T() const {
        return val;
    }

    TrackedVar(const T &v) {
        val = v;
        values.insert(val);
    }

#define DEF_BOP_THAT_ASSIGNS(OP)                                                                   \
    TrackedVar &operator OP(const T & rhs) {                                                       \
        val OP rhs;                                                                                \
        values.insert(val);                                                                        \
        return *this;                                                                              \
    }                                                                                              \
    TrackedVar &operator OP(const TrackedVar & rhs) {                                              \
        val OP rhs.val;                                                                            \
        values.insert(val);                                                                        \
        return *this;                                                                              \
    }
    DEF_BOP_THAT_ASSIGNS(=)
    DEF_BOP_THAT_ASSIGNS(*=)
    DEF_BOP_THAT_ASSIGNS(/=)
    DEF_BOP_THAT_ASSIGNS(%=)
    DEF_BOP_THAT_ASSIGNS(+=)
    DEF_BOP_THAT_ASSIGNS(-=)
    DEF_BOP_THAT_ASSIGNS(<<=)
    DEF_BOP_THAT_ASSIGNS(>>=)
    DEF_BOP_THAT_ASSIGNS(&=)
    DEF_BOP_THAT_ASSIGNS(^=)
    DEF_BOP_THAT_ASSIGNS(|=)

#undef DEF_BOP_THAT_ASSIGNS

    // pre inc/dec, return a reference like for built-in types: `T &r = ++x;`
    TrackedVar &operator++() {
        ++val;
        values.insert(val);
        return *this;
    }
    TrackedVar &operator--() {
        --val;
        values.insert(val);
        return *this;
    }
    // post inc/dec
    T operator++(int) {
        T ret = val++;
        values.insert(val);
        return ret;
    }
    T operator--(int) {
        T ret = val--;
        values.insert(val);
        return ret;
    }
};