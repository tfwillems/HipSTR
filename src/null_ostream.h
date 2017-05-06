#ifndef NULL_OSTREAM_H_
#define NULL_OSTREAM_H_

#include <ostream>

class NullOstream : public std::ostream {
 private:

  class NullOutputBuffer : public std::streambuf {
  public:
    virtual int overflow(int c = EOF){
      return 1;
    }

    virtual std::streamsize xsputn(const char* s, std::streamsize n){
      return n;
    }
  };

  NullOutputBuffer out_;

 public:
 NullOstream() : std::ostream(&out_){}
};

#endif
