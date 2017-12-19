#include "vcf_writer.h"

bool tuple_comparator(RecordTuple* r1, RecordTuple* r2){
  return r1->pos() > r2->pos();
}

void VCFWriter::add_vcf_record(const std::string& chrom, int32_t record_pos, const std::string& record_text){
  if (!open_)
    printErrorAndDie("Cannot invoke add_vcf_record() on a non-open VCFWriter");

  // If we're changing chromosomes, output all the current records
  if (chrom.compare(chrom_) != 0){
    write_all_records();
    chrom_ = chrom;
  }
  else {
    // Output any existing records that definitely precede all future records
    while (!record_heap_.empty()){
      std::pop_heap(record_heap_.begin(), record_heap_.end(), tuple_comparator);
      RecordTuple* best = record_heap_.back(); record_heap_.pop_back();
      if (best->pos() < record_pos - MAX_RECORD_PAD){
	str_vcf_ << best->text() << std::endl;
	delete best;
      }
      else {
	record_heap_.push_back(best);
	std::push_heap(record_heap_.begin(), record_heap_.end(), tuple_comparator);
	break;
      }
    }
  }

  // Add the newest record to the heap
  record_heap_.push_back(new RecordTuple(record_pos, record_text));
  std::push_heap(record_heap_.begin(), record_heap_.end(), tuple_comparator);
}
