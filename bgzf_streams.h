#ifndef __BGZF_STREAMS_H__
#define __BGZF_STREAMS_H__

/***********************************************************************

Copyright (C) 2015 Assaf Gordon (agordon@nygenome.org)
License: GPLv2-or-later.

STl/ostream and STL/istream wrappers around HTSLIB's BGZF module.

Typical usage:

	bgzfostream a("1.bgzf");

	a << "hello world";
	a << 42;
	a << "\nanother test" << endl;


bgzfostream inherits from std::ostream and can be used interchangeably.
bgzfistream inherits from std::istream and can be used interchangeably.

TODO:
`. Replace 'err()' with proper STL exceptions.

************************************************************************/

#include <err.h>
#include <ostream>
#include <stdexcept>

#include "vcflib/tabixpp/htslib/htslib/bgzf.h"

class bgzf_streambuf : public std::streambuf {
 private:
  BGZF* _fp;
  std::string filename;
  int cur_val;

 public:
 bgzf_streambuf(): _fp(NULL){ 
    cur_val = -999;
  }
  
  virtual ~bgzf_streambuf(){
    close();
  }
  
  void open(const char *_filename, const char *mode){
    if (_fp != NULL)
      throw std::invalid_argument("bgzf_streambuf: open: called on an open stream");
    
    _fp = bgzf_open(_filename, mode);
    if (_fp == NULL)
      err(1,"bgzf_open(%s,%s) failed", _filename, mode);
    filename = _filename;
  }
  
  void close(){
    if (_fp == NULL)
      return;
    
    int i = bgzf_close(_fp);
    if (i != 0)
      err(1,"bgzf_close(%s) failed", filename.c_str());
    
    _fp = NULL;
    filename = "";
  }
  
  virtual int uflow(){
    if (cur_val != -999){
      int res = cur_val;
      cur_val = -999;
      return res;
    }
    int res = bgzf_getc(_fp);
    if (res == -2) err(1, "bgzf_getc() failed");
    if (res == -1) return EOF;
    return res;
  }
  
  virtual int underflow(){
    if ( _fp == NULL)
      throw std::invalid_argument("bgzf_streambuf: underflow: called on non-open stream");
    if (cur_val != -999)
      return cur_val;
    cur_val = bgzf_getc(_fp);
    if (cur_val == -2) err(1, "bgzf_getc() failed");
    if (cur_val == -1) cur_val = EOF;
    return cur_val;
  }

  virtual int overflow(int c = EOF){
    char z = (char) c;
    ssize_t i = bgzf_write(_fp, &z, 1);
    if (i < 0) {
      err(1,"bgzf_write(%s) failed", filename.c_str());
      return EOF;
    }
    return c;
  }

  virtual std::streamsize xsputn (const char* s, std::streamsize n){
    if ( _fp == NULL)
      throw std::invalid_argument("bgzf_streambuf: overflow: called on non-open stream");
    
    ssize_t i = bgzf_write(_fp, s, n);
    if (i < 0) {
      err(1,"bgzf_write(%s) failed", filename.c_str());
      return EOF;
    }
    if ((std::streamsize)i != n) {
      err(1,"bgzf_write(%s) wrote only %zd, asked for %zu bytes",
	  filename.c_str(), i, n);
      return EOF;
    }
    return (std::streamsize)n;
  }
};

class bgzfistream : public std::istream {
 protected:
  bgzf_streambuf buf;
 public:
  
 bgzfistream(const char* filename, const char *mode="r") : std::istream(0) {
    buf.open(filename, mode);
    rdbuf(&buf);
  }
  
 bgzfistream() : std::istream(0) {}
  
  void open(const char* filename, const char *mode="r") {
    buf.open(filename, mode);
    rdbuf(&buf);
  }

  void close(){
    buf.close();
  }
};

class bgzfostream : public std::ostream {
 protected:
  bgzf_streambuf buf;
 public:
  
 bgzfostream(const char* filename, const char *mode="w") : std::ostream(0) {
    buf.open(filename, mode);
    rdbuf(&buf);
  }
  
 bgzfostream() : std::ostream(0) {}
  
  void open(const char* filename, const char *mode="w") {
    buf.open(filename, mode);
    rdbuf(&buf);
  }

  void close(){
    buf.close();
  }
};

#endif
