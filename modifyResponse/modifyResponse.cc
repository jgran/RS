#include "modifyResponse.h"
#include "do_convolution.h"
#include "do_scale.h"

int main(int argc, char **argv) {

  std::vector<std::string> allowed_modMethod;
  allowed_modMethod.push_back("scale");
  allowed_modMethod.push_back("convolution");

  std::vector<std::string> allowed_modType;
  allowed_modType.push_back("low_tail");
  allowed_modType.push_back("high_tail");
  allowed_modType.push_back("low_and_high_tail");
  allowed_modType.push_back("core");

  if (argc != 4) {
    std::cout << "USAGE: modifyResponse <modMethod> <modType> <modValue>" << std::endl;
    std::cout << std::endl;
    std::cout << "Allowed modMethods:" << std::endl;
    for(std::vector<std::string>::const_iterator it = allowed_modMethod.begin(); it != allowed_modMethod.end(); it++) std::cout << *it << std::endl;
    std::cout << std::endl;
    std::cout << "Allowed modTypes:" << std::endl;
    for(std::vector<std::string>::const_iterator it = allowed_modType.begin(); it != allowed_modType.end(); it++) std::cout << *it << std::endl;
    std::cout << std::endl;
    return 1;
  }

  std::string modMethod = argv[1];
  std::string modType = argv[2];
  float modValue = std::atof(argv[3]);
  
  bool good_modMethod = false;
  for(std::vector<std::string>::const_iterator it = allowed_modMethod.begin(); it != allowed_modMethod.end(); it++){
    if(modMethod == *it){
      good_modMethod = true;
      break;
    } 
  }
  if(! good_modMethod){
    std::cout << "modMethod: " << modMethod << " is not supported." << std::endl;
    return 1;
  }

  bool good_modType = false;
  for(std::vector<std::string>::const_iterator it = allowed_modType.begin(); it != allowed_modType.end(); it++){
    if(modType == *it){
      good_modType = true;
      break;
    } 
  }
  if(! good_modType){
    std::cout << "modType: " << modType << " is not supported." << std::endl;
    return 1;
  }

   
  if(modMethod == "convolution"){
    convol myconvol;
    myconvol.do_convol(modType, modValue);
  }
  else if(modMethod == "scale"){
    scale myscale;
    myscale.do_scale(modType, modValue);
  }
  else{
    std::cout << "shouldn't get here!" << std::endl;
  }

  return 0;
}
