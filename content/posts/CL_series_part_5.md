---
title: "Essential CLI Tools You Should Know - Part 5: bottom & bat"
author: "Badran Elshenawy"
date: 2025-01-25
categories:
  - "CLI Tools"
  - "Linux"
  - "Productivity"
tags:
  - "CLI"
  - "Linux"
  - "SysAdmin"
  - "DevOps"
  - "Automation"
description: "A deep dive into two essential CLI tools, `bottom` and `bat`, that will enhance your workflow."
slug: "cli-tools-bottom-bat"
draft: false
output: hugodown::md_document
aliases:
  - "/posts/CL_series_part_5/"
summary: "Explore `bottom`, a better `htop`, and `bat`, a modern `cat` replacement to improve your Linux workflow."
featured: true
rmd_hash: ea8169414cbca5ef

---

# 🛠️ Essential CLI Tools You Should Know - Part 5: `bottom` & `bat`

Monitoring system resources and viewing files are **daily tasks** in the command line. But default tools like `htop` and `cat` **lack modern usability features**. Thankfully, two **powerful replacements** make these tasks **faster, clearer, and more efficient**:

-   **`bottom`** → A feature-rich `htop` alternative for system monitoring

-   **`bat`** → A smarter, more user-friendly replacement for `cat`

------------------------------------------------------------------------

## 🌟 `bottom` \-- A Modern `htop` Alternative

If you're still using `htop` for system monitoring, **you're missing out**. `bottom` (`btm`) is a **lightweight, high-performance** process monitor with real-time graphs and a beautiful interface.

### 🔹 Why `bottom`?

✅ **Faster than `htop`** -- Rust-powered for efficiency  
✅ **Beautiful UI** -- Real-time CPU, memory, disk, and network graphs  
✅ **Tree view for processes** -- See parent-child relationships  
✅ **Customizable widgets** -- Tailor your monitoring experience  
✅ **Battery monitoring** -- Useful for laptop users

### ⚡ Quick Start

#### 📥 Install:

-   **Ubuntu/Debian:** `sudo apt install bottom`

-   **MacOS:** `brew install bottom`

-   **Rust Users:** `cargo install bottom`

#### 💡 Common Usage:

    # Start system monitor
    btm  

    # Minimal UI with essential stats
    btm --basic  

    # View processes in hierarchical tree view
    btm --tree  

🔥 **Why switch?** `bottom` gives a **better, more insightful view** of system resources with **less overhead** than `htop`.

------------------------------------------------------------------------

## 🌟 `bat` \-- The Modern `cat` Replacement

The default `cat` command works, but it **lacks modern usability features**. `bat` makes file viewing **more powerful and readable** with syntax highlighting, line numbers, and Git integration.

### 🔹 Why `bat`?

✅ **Syntax highlighting** -- Works for code, configs, logs  
✅ **Git integration** -- Shows diffs and untracked changes  
✅ **Built-in paging** -- Scroll through large files automatically  
✅ **Supports multiple files** -- Compare files side by side  
✅ **Drop-in replacement for `cat`** -- Works the same way, just better

### ⚡ Quick Start

#### 📥 Install:

-   **Ubuntu/Debian:** `sudo apt install bat`

-   **MacOS:** `brew install bat`

#### 💡 Common Usage:

    # View file with syntax highlighting
    bat filename  

    # Disable extra formatting for raw output
    bat --style=plain filename  

    # Compare multiple files side by side
    bat file1 file2  

    # Use bat in pipelines
    cat file | bat  

🔥 **Why switch?** `bat` makes reading files **visually appealing** and **more functional** while keeping `cat`'s simplicity.

------------------------------------------------------------------------

## ✅ Final Thoughts

Both `bottom` and `bat` **replace outdated CLI tools with modern, feature-rich alternatives** that improve productivity.

### 🔥 Pro Tips:

🔸 **Use `btm --tree`** to quickly analyze resource-hungry processes  
🔸 **Set `bat` as your default pager** by adding this to your shell config:

    export PAGER="bat"

🔸 **Combine `fd` + `bat`** for rapid file searches:

    fd "pattern" | bat --paging=never

🚀 **These tools will supercharge your workflow!**

