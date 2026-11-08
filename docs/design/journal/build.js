// Inline the shared CSS/JS into each screen so every HTML file stands alone.
const fs = require('fs'), path = require('path');
const src = path.join(__dirname, 'src');
const css = fs.readFileSync(path.join(src, '_journal.css'), 'utf8');
const js = fs.readFileSync(path.join(src, '_journal.js'), 'utf8');
for (const f of fs.readdirSync(src).filter(f => /^\d\d-.*\.html$/.test(f))) {
  const h = fs.readFileSync(path.join(src, f), 'utf8')
    .replace('<!--@BASE-->', `<style>\n${css}\n</style>\n<script>\n${js}\n</script>`);
  fs.writeFileSync(path.join(__dirname, 'html', f), h);
  console.log('built', f);
}
