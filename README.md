## To test

0. Misc updates

Update defaults in public/`manifest.json`
Change favicons and logos as needed

1. Install dependencies (done once)

```
npm install
```

Update `package.json`

2. Run app to test locally

Runs the app in the development mode.
Open http://localhost:3000 to view it in your browser.

```
npm start
```

Another host is needed to serve files from the `config.json` base url. For example, if in R:
```
base_url <- "http://localhost:8000/"
```

Then a localhost server needs to be launched:
```
http-server ./ --cors -p 8000
```

Alternatively push files to AWS or Github and reference them directly:
```
base_url <- "https://raw.githubusercontent.com/JEFworks/vitessce..."
```

3. Build and rename

```
npm run build
mv build/ docs/
```

Builds the app for production to the build folder.
Rename to docs for Github pages.