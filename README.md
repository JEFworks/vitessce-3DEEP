## To test

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

Note, due to caching, seems like another host is needed to serve files:
```
http-server ./ --cors -p 8000
```

3. Build and rename

```
npm run build
mv build/ docs/
```

Builds the app for production to the build folder.
Rename to docs for Github pages.