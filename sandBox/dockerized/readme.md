From the directory containing /app, Dockerfile, config file and the bash script, run <br><br>
```docker build -t tagname .```

Once the image is built, initiate the app using:<br><br>
```docker run -p 3838:80 tagname```

