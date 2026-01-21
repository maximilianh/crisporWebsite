int read_and_increment_file_counter(const char*const filename, int startValue=0, int maxValue=INT32_MAX) {
  ifstream in(filename);
  int value;
  if (!(in >> value) || value > maxValue) {
    value = startValue;
  }
  ofstream out(filename, ios_base::out | ios_base::trunc);
  out << (value+1) << "\n";
  return value;
}

static std::string af_counter = std::to_string(read_and_increment_file_counter("accessfold.counter.tmp",1, 7));

string getDebugFileName(const string& basename) {
  auto nanos = std::chrono::system_clock::now().time_since_epoch();
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(nanos);
  string name = basename; name += '_';
  const char *envp = std::getenv("DEBUG_ID");
  if (envp!=NULL) {
    name+=envp; name+='_';
  }
  name.append(af_counter).append(1, '_');
  name += std::to_string(ms.count());
  name += ".txt";
  return name;
}