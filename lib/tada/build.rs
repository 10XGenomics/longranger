use std::process::Command;
use std::env;

fn main() {

    let _path = env::var_os("OUT_DIR").unwrap();
    let path = _path.to_str().unwrap();

    let cmd = format!("git describe --tags --dirty > {}/version.txt", path);
    let output = Command::new("sh")
                        .arg("-c")
                        .arg(cmd)
                        .output()
                        .expect("failed to execute process");
}