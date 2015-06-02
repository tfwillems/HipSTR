#ifndef __BGZF_OSTREAM_H__
#define __BGZF_OSTREAM_H__

/***********************************************************************

Copyright (C) 2015 Assaf Gordon (agordon@nygenome.org)
License: GPLv2-or-later.

An STl/Ostream wrapper around HTSLIB's BGZF module.

Typical usage:

	bgzfostream a("1.bgzf");

	a << "hello world";
	a << 42;
	a << "\nanother test" << endl;


bgzfostream inherits from 'std::ostream' and can be used interchangibly.


TODO:
1. implement read functionality as well, rename 'bgzf_write_streambuf'
   to 'bgzf_streambuf'.
2. replace 'err()' with proper STL exceptions.

************************************************************************/

#include <err.h>
#include <ostream>
#include <stdexcept>

#include "vcflib/tabixpp/htslib/htslib/bgzf.h"

class bgzf_write_streambuf : public std::streambuf {
 private:
  BGZF* _fp;
  std::string filename;
  
 public:
 bgzf_write_streambuf(): _fp(NULL){ }
  
  virtual ~bgzf_write_streambuf(){
    close();
  }
  
  void open(const char *_filename, const char *mode){
    if (_fp != NULL)
      throw std::invalid_argument("bgzf_write_streambuf: open: called on an open stream");
    
    _fp = bgzf_open ( _filename, mode );
    if (_fp == NULL )
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
  
  virtual int overflow(int c = EOF){
    if ( _fp == NULL)
      throw std::invalid_argument("bgzf_write_streambuf: overflow: called on non-open stream");
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
      throw std::invalid_argument("bgzf_write_streambuf: overflow: called on non-open stream");
    
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

class bgzfostream : public std::ostream {
 protected:
  bgzf_write_streambuf buf;
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
