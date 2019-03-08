#![feature(proc_macro_hygiene, decl_macro)]

#[macro_use] extern crate rocket;

#[cfg(test)] mod tests;



#[get("/lisbert/<verb>")]
fn hello(verb: String) -> String {format!("Lisa {}!", verb)}

fn main() {
    rocket::ignite().mount("/", routes![hello]).launch();
}
