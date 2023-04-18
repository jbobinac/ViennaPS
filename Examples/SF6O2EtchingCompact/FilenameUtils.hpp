#pragma once

#include <regex>
#include <string>
#include <utility>

std::pair<std::string, std::string>
//std::string
extractParameters(const std::string &filename,
                  const std::string &regexPattern) {
  std::regex rgx(regexPattern);
  std::smatch smatch;

  if (std::regex_search(filename.begin(), filename.end(), smatch, rgx)) {
    if (smatch.size() < 2)
      return {};
    return {smatch[1], smatch[2]};
  } else {
    return {};
  }
}
