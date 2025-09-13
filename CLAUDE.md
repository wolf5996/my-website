# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a personal website for Badran Elshenawy, a bioinformatician and postdoc researcher. The site is built with Hugo using the Coder theme and deployed on Netlify. It features a personal blog with 128+ posts covering bioinformatics, R programming, command-line tools, and computational biology topics.

## Development Commands

### Local Development
```bash
# Start local development server with drafts
hugo server -D

# Start server without drafts
hugo server

# Build for production
hugo
```

### Content Management
```bash
# Create new blog post
hugo new posts/your-post-title.md

# Create new content with archetype
hugo new content/about/new-page.md
```

## Site Configuration

- **Main config**: `hugo.toml` - Contains site metadata, theme configuration, social links, and navigation
- **Build config**: `netlify.toml` - Netlify deployment settings with Hugo v0.141.0
- **Theme**: Uses the "Coder" theme located in `themes/Coder/`

## Content Structure

- `content/posts/` - Blog posts (128+ markdown files covering bioinformatics tutorials and series)
- `content/about.md` - Personal biography and expertise overview
- `static/` - Static assets (images, etc.)
- `public/` - Generated site output (Hugo build target)

## Key Features

- Personal blog focusing on bioinformatics, R programming, and computational biology
- Multi-series content including DESeq2, tidyverse, foundational genomics, and bulk RNA-seq tutorials
- Social media integration (LinkedIn, GitHub, Oxford profile, Twitter)
- Light/dark mode support via Coder theme
- Syntax highlighting enabled for code blocks
- RSS feed capabilities

## Deployment

The site is automatically deployed to Netlify when changes are pushed to the main branch. The live site is available at: https://badran-elshenawy.netlify.app/

## Content Guidelines

Blog posts are written in Markdown and typically cover:
- Bioinformatics workflows and tutorials
- R programming and tidyverse techniques  
- Command-line tools and automation
- Git and reproducibility practices
- High-performance computing for biology

When creating new content, follow the existing naming conventions in `content/posts/` and include appropriate front matter for Hugo processing.