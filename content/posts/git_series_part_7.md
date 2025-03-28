---
title: "Version Control Series – Part 7: Git Diff & Log – Tracking Changes Like a Pro"
author: "Badran Elshenawy"
date: 2025-02-04T12:30:00Z
categories:
  - "Version Control"
  - "Git"
  - "Bioinformatics"
tags:
  - "Git"
  - "Version Control"
  - "Bioinformatics"
  - "Collaboration"
  - "Reproducibility"
description: "A comprehensive guide on using git diff and git log to track changes in your repository, with practical examples and best practices."
slug: "git-diff-log-tracking"
draft: false
output: hugodown::md_document
aliases:
  - "/posts/git_diff_log_tracking/"
summary: "Master git diff and git log to track changes efficiently and enhance your version control workflow."
featured: true
rmd_hash: 1f02a666effbb1b1

---

# 🪠 Version Control Series -- Post 7: Git Diff & Log -- Tracking Changes Like a Pro

## 🔍 Why Tracking Changes Matters

Version control isn't just about saving files---it's about knowing **what changed, when, and why**. Git provides powerful tools to track modifications, helping you debug issues, review work, and collaborate effectively.

One of the biggest challenges in software development and data science is **tracking changes across different versions of your code**. With Git, you can pinpoint exactly what has changed between commits, branches, and even contributors. In this guide, we'll explore how to track changes efficiently using `git diff` and `git log` and how modern tools like GitLens and Git Graph simplify the process.

------------------------------------------------------------------------

## 🔀 1️⃣ Viewing Changes with `git diff`

The `git diff` command allows you to see **line-by-line differences** between file versions. This is useful before committing changes or reviewing differences between branches.

### ✅ Basic `git diff` Usage

``` bash
git diff
```

This command shows **unstaged changes**, highlighting modified lines since the last commit.

### ✅ Viewing Staged Changes

``` bash
git diff --staged
```

Use this to see changes that are staged for the next commit.

### ✅ Comparing with the Latest Commit

``` bash
git diff HEAD
```

Compares the working directory with the latest commit, helping you track what has changed since your last save point.

### ✅ Comparing Two Commits

``` bash
git diff commit1 commit2
```

Replace `commit1` and `commit2` with actual commit hashes to see differences between any two points in history.

> 🔹 **Pro Tip:** Use `git diff --color-words` for a more readable diff that highlights word changes instead of whole lines.

------------------------------------------------------------------------

## 📜 2️⃣ Checking Commit History with `git log`

The `git log` command helps you explore **the history of your project**.

### ✅ Viewing Basic Commit History

``` bash
git log
```

Lists commits in reverse chronological order.

### ✅ Simplified Log for Quick Overview

``` bash
git log --oneline --graph --all
```

This provides a **clean, visual representation** of the commit tree, making it easier to understand the history of your repository.

### ✅ Viewing the Last N Commits with Diffs

``` bash
git log -p -2
```

Displays the **last two commits** with detailed changes.

> 🔹 **Pro Tip:** Customize your log view with `git log --pretty=format:'%h - %an, %ar : %s'` for a cleaner history output.

------------------------------------------------------------------------

## 📊 3️⃣ Making Git History Easier with GitLens & Git Graph

While Git's command-line tools are powerful, **VS Code extensions like GitLens and Git Graph** make tracking changes much easier.

### 🔹 **GitLens** -- Enhanced Git Capabilities

[GitLens](https://marketplace.visualstudio.com/items?itemName=eamodio.gitlens) integrates Git deeply into VS Code, allowing you to: - View commit details **inline** - See who last modified a line of code - Navigate commit history effortlessly

### 🔹 **Git Graph** -- A Visual Git History

[Git Graph](https://marketplace.visualstudio.com/items?itemName=mhutchie.git-graph) offers: - A **graphical representation** of branches and commits - The ability to **easily switch between branches** - Quick visualization of **merge histories**

> 🔹 **Pro Tip:** Combine GitLens with Git Graph in VS Code for the ultimate Git experience!

------------------------------------------------------------------------

## 🎯 Key Takeaways

✅ **Use `git diff` to review changes before committing.**  
✅ **Use `git log` to navigate commit history effectively.**  
✅ **Extensions like GitLens & Git Graph simplify Git tracking in VS Code.**

📌 **Next up: Best Practices for Writing Commit Messages!** Stay tuned! 🚀

👇 How do you track changes in your Git workflow? Let's discuss!  
#Git #VersionControl #Bioinformatics #Reproducibility #OpenScience #ComputationalBiology

