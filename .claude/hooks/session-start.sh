#!/bin/bash
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

install_plugin() {
  local marketplace="$1" plugin="$2"
  if ! claude plugin list 2>/dev/null | grep -q "$plugin"; then
    claude plugin marketplace add "$marketplace" || return 0
    claude plugin install "$plugin" || true
  fi
}

install_plugin ayghri/i-have-adhd i-have-adhd@i-have-adhd
install_plugin adaptyvbio/protein-design-skills adaptyv@protein-design-skills
install_plugin mims-harvard/ToolUniverse tooluniverse@tooluniverse
install_plugin AlterLab-IEU/AlterLab-Academic-Skills alterlab-writing-tools@alterlab-academic-skills
install_plugin AlterLab-IEU/AlterLab-Academic-Skills alterlab-data-science@alterlab-academic-skills
install_plugin orchestra-research/AI-research-SKILLs fine-tuning@ai-research-skills
install_plugin orchestra-research/AI-research-SKILLs evaluation@ai-research-skills
install_plugin orchestra-research/AI-research-SKILLs inference-serving@ai-research-skills
install_plugin orchestra-research/AI-research-SKILLs prompt-engineering@ai-research-skills
install_plugin orchestra-research/AI-research-SKILLs agents@ai-research-skills

CONFIG_DIR="${CLAUDE_CONFIG_DIR:-$HOME/.claude}"
mkdir -p "$CONFIG_DIR"
touch "$CONFIG_DIR/.i-have-adhd-always"
