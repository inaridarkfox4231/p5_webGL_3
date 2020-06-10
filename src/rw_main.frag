void mainImage(out vec4 fragColor, in vec2 fragCoord){
  vec2 p = fragCoord.xy / iResolution.xy;

  // Time varying pixel color
  vec4 col = texture(iChannel0, p);
  vec4 rw = texture(iChannel1, p);
  col = mix(col, rw, rw.a);

  // Output to screen
  fragColor = col;
}
