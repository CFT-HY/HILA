
template <typename T>
class tst {
public:
  field<T> k;
  void sub(field<T> &a, const field<T> &b, parity p)
  {
    a[p] -= b[X];
  }
};
