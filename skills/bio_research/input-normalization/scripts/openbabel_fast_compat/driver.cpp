// SPDX-License-Identifier: GPL-2.0-only
// Narrow CLI adapter for the pinned Open Babel operation in gen3d_fast.cpp.
#include <openbabel/obconversion.h>
#include <openbabel/op.h>

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  std::string input_format, input_path, output_path;
  bool sdf_output = false, hydrogens = false, fast = false;
  for (int i = 1; i < argc; ++i) {
    const std::string argument(argv[i]);
    if (argument.size() > 2 && argument.compare(0, 2, "-i") == 0) {
      if (!input_format.empty() || i + 1 >= argc) return 2;
      input_format = argument.substr(2);
      input_path = argv[++i];
    } else if (argument == "-osdf") {
      sdf_output = true;
    } else if (argument == "-O" && i + 1 < argc) {
      output_path = argv[++i];
    } else if (argument == "-h") {
      hydrogens = true;
    } else if (argument == "--gen3d" && i + 1 < argc) {
      fast = std::string(argv[++i]) == "fast";
    } else {
      std::cerr << "Unsupported fast compatibility argument: " << argument << '\n';
      return 2;
    }
  }
  if (input_format.empty() || input_path.empty() || output_path.empty() ||
      !sdf_output || !hydrogens || !fast) {
    std::cerr << "Expected -iFORMAT INPUT -osdf -O OUTPUT -h --gen3d fast\n";
    return 2;
  }

  std::ifstream input(input_path.c_str(), std::ios::binary);
  std::ofstream output(output_path.c_str(), std::ios::binary);
  if (!input || !output) {
    std::cerr << "Cannot open the conversion input or output\n";
    return 1;
  }
  OpenBabel::OBConversion conversion(&input, &output);
  if (!conversion.SetInAndOutFormats(input_format.c_str(), "sdf") ||
      !OpenBabel::OBOp::FindType("sbddgen3dfast")) {
    std::cerr << "Open Babel format or fast compatibility operation unavailable\n";
    return 1;
  }
  // Keep Open Babel's native conversion/transform pipeline, including operation
  // ordering, hydrogen addition, record titles and multi-record error behavior.
  conversion.AddOption("h", OpenBabel::OBConversion::GENOPTIONS);
  conversion.AddOption("sbddgen3dfast", OpenBabel::OBConversion::GENOPTIONS, "fast");
  const int converted = conversion.Convert();
  output.flush();
  if (!output) {
    std::cerr << "Failed writing converted SDF\n";
    return 1;
  }
  std::cerr << converted << (converted == 1 ? " molecule converted\n" : " molecules converted\n");
  return 0;
}
