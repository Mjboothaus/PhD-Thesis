import markdown
from markdown_include.include import MarkdownInclude

# Markdown Extensions

markdown_include = MarkdownInclude(
    configs={'base_path': '/srv/content/', 'encoding': 'iso-8859-1'}
)
html = markdown.markdown(source, extensions=[markdown_include])
