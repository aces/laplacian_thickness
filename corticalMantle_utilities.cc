#include "corticalMantle.h"

void corticalMantle::outputVolume(STRING filename) {
  save_label_volume(filename, this->clsFile, this->mantle, 255);
}

void corticalMantle::outputVolume() {
  save_label_volume(this->outputFile, this->clsFile, this->mantle, 255);
}
